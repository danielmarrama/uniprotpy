#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')

import re
import os
import pandas as pd
import requests


class ProteomeSelector:
  def __init__(self, taxon_id):
    self.taxon_id = taxon_id
    
    # get proteome list for species and count number of proteomes
    self.proteome_list = self._get_proteome_list()
    self.num_of_proteomes = len(self.proteome_list) + 1 # +1 because "all proteins" is also a candidate proteome

  def select_proteome(self):
    """
    Select the best proteome to use for a species. Return the proteome ID, 
    proteome taxon, and proteome type.
    
    Check UniProt for all candidate proteomes associated with that
    taxon. Then do the following checks:

    1. Are there any representative proteomes?
    2. Are there any reference proteomes?
    3. Are there any non-redudant proteomes?
    4. Are there any other proteomes?

    If yes to any of the above, check if there are ties. If there are ties,
    then select the proteome with the most proteins.

    If no to all of the above, then get every protein associated with
    the taxon ID using the get_all_proteins method.
    """
    # if species_dir already exists then return the already selected proteome, else create dir
    if os.path.exists(f'./data/{self.taxon_id}'):
      print(f'Proteome already selected for {self.taxon_id}.')
      return []
    else:
      os.makedirs(f'./data/{self.taxon_id}')

    # if there is no proteome_list, get all proteins associated with that taxon ID
    if self.proteome_list.empty:
      self._get_all_proteins()
      return 'None', self.taxon_id, 'All-proteins'

    if self.proteome_list['isRepresentativeProteome'].any():
      proteome_type = 'Representative'
      self.proteome_list = self.proteome_list[self.proteome_list['isRepresentativeProteome']]
      proteome_id, proteome_taxon = self._get_proteome_with_most_proteins()
      # self._get_gp_proteome_to_fasta(proteome_id, proteome_taxon)
    
    elif self.proteome_list['isReferenceProteome'].any():
      proteome_type = 'Reference'
      self.proteome_list = self.proteome_list[self.proteome_list['isReferenceProteome']]
      proteome_id, proteome_taxon = self._get_proteome_with_most_proteins()
      # self._get_gp_proteome_to_fasta(proteome_id, proteome_taxon)

    elif 'redundantTo' not in self.proteome_list.columns:
      proteome_type = 'Other'
      proteome_id, proteome_taxon = self._get_proteome_with_most_proteins()
    
    elif self.proteome_list['redundantTo'].isna().any():
      proteome_type = 'Non-redundant'
      self.proteome_list = self.proteome_list[self.proteome_list['redundantTo'].isna()]
      proteome_id, proteome_taxon = self._get_proteome_with_most_proteins()
    
    else:
      proteome_type = 'Other'
      proteome_id, proteome_taxon = self._get_proteome_with_most_proteins()
    
    # sanity check to make sure proteome.fasta is not empty
    if os.stat(f'./data/{self.taxon_id}/proteome.fasta').st_size == 0:
      proteome_id = 'None'
      proteome_taxon = self.taxon_id
      proteome_type = 'All-proteins'
      self._get_all_proteins()

    proteome_data = [proteome_id, proteome_taxon, proteome_type]
    return proteome_data

  def proteome_to_csv(self):
    """
    Write the proteome data for a species to a CSV file for later use.
    """
    from Bio import SeqIO

    # read in the FASTA file and then get the gene priority IDs if they exist
    proteins = list(SeqIO.parse(f'data/{self.taxon_id}/proteome.fasta', 'fasta'))

    gp_proteome_path = f'data/{self.taxon_id}/gp_proteome.fasta'
    if os.path.isfile(gp_proteome_path):
      gp_ids = [str(protein.id.split('|')[1]) for protein in list(SeqIO.parse(gp_proteome_path, 'fasta'))]
    else:
      gp_ids = []

    # start collecting proteome data
    proteome_data = []
    for protein in proteins:
      uniprot_id = protein.id.split('|')[1]
      gp = 1 if uniprot_id in gp_ids else 0

      # TODO: look into using HUGO to map old gene names to new ones

      try:
        gene = re.search('GN=(.*?) ', protein.description).group(1)
      except AttributeError:
        try: # some gene names are at the end of the header
          gene = re.search('GN=(.*?)$', protein.description).group(1)
        except AttributeError:
          gene = ''
      try: # protein existence level
        pe_level = int(re.search('PE=(.*?) ', protein.description).group(1))
      except AttributeError:
        pe_level = 0
      
      proteome_data.append([protein.id.split('|')[0], gene, uniprot_id, gp, pe_level, str(protein.seq)])
    
    columns = ['Database', 'Gene Symbol', 'UniProt ID', 'Gene Priority', 'Protein Existence Level', 'Sequence']
    pd.DataFrame(proteome_data, columns=columns).to_csv(f'data/{self.taxon_id}/proteome.csv', index=False)

  def _get_proteome_list(self):
    """
    Get a list of proteomes for a species from the UniProt API.
    Check for proteome_type:1 first, which are the representative or
    reference proteomes.

    If there are no proteomes, return empty DataFrame.
    """
    # URL to get proteome list for a species - use proteome_type:1 first
    url = f'https://rest.uniprot.org/proteomes/stream?format=xml&query=(proteome_type:1)AND(taxonomy_id:{self.taxon_id})'
    
    try:
      proteome_list = pd.read_xml(requests.get(url).text)
    except ValueError:
      try: # delete proteome_type:1 from URL and try again
        url = url.replace('(proteome_type:1)AND', '')
        proteome_list = pd.read_xml(requests.get(url).text)
      except ValueError: # if there are no proteomes, return empty DataFrame
        return pd.DataFrame()

    # remove the namespace from the columns
    proteome_list.columns = [x.replace('{http://uniprot.org/proteome}', '') for x in proteome_list.columns]
    return proteome_list

  def _get_all_proteins(self):
    """
    Get every protein associated with a taxon ID on UniProt.
    Species on UniProt will have a proteome, but not every protein is
    stored within those proteomes. There is a way to get every protein
    using the taxonomy part of UniProt. 
    """
    # URL link to all proteins for a species - size = 500 proteins at a time
    url = f'https://rest.uniprot.org/uniprotkb/search?format=fasta&'\
          f'query=taxonomy_id:{self.taxon_id}&size=500' 

    # loop through all protein batches and write proteins to FASTA file
    for batch in self._get_protein_batches(url):
      with open(f'data/{self.taxon_id}/proteome.fasta', 'a') as f:
        f.write(batch.text)

  def _get_protein_batches(self, batch_url):
    """
    Get a batch of proteins from UniProt API because it limits the
    number of proteins you can get at once. Yield each batch until the 
    URL link is empty.
    
    Args:
      batch_url (str): URL to get all proteins for a species.
    """
    while batch_url:
      r = requests.get(batch_url)
      r.raise_for_status()
      yield r
      batch_url = self._get_next_link(r.headers)

  def _get_next_link(self, headers):
    """
    UniProt will provide a link to the next batch of proteins when getting
    all proteins for a species' taxon ID.
    We can use a regular expression to extract the URL from the header.

    Args:
      headers (dict): Headers from UniProt API response.
    """
    re_next_link = re.compile(r'<(.+)>; rel="next"') # regex to extract URL
    if 'Link' in headers:
      match = re_next_link.match(headers['Link'])
      if match:
        return match.group(1)

  def _get_gp_proteome_to_fasta(self, proteome_id, proteome_taxon):
    """
    Write the gene priority proteome to a file. 
    This is only for representative and reference proteomes.
    Depending on the species group, the FTP URL will be different.

    Args:
      proteome_id (str): Proteome ID.
      proteome_taxon (str): Taxon ID for the proteome.
    """
    import gzip
    
    group = self.species_df[self.species_df['Taxon ID'].astype(str) == self.taxon_id]['Group'].iloc[0]
    ftp_url = f'https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/'
    
    if group == 'archeobacterium':
      ftp_url += f'Archaea/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif group == 'bacterium':
      ftp_url += f'Bacteria/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif group in ['vertebrate', 'other-eukaryote']:
      ftp_url += f'Eukaryota/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif group in ['virus', 'small-virus', 'large-virus']:
      ftp_url += f'Viruses/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    
    r = requests.get(ftp_url, stream=True)
    try:
      r.raise_for_status()
    except:
      return

    # unzip the request and write the gene priority proteome to a file
    with open(f'data/{self.taxon_id}/gp_proteome.fasta', 'wb') as f:
      f.write(gzip.open(r.raw, 'rb').read())

  def _get_proteome_to_fasta(self, proteome_id):
    """
    Get the FASTA file for a proteome from UniProt API.
    Include all isoforms and do not compress the file.
    """
    url = f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&includeIsoform=true&query=(proteome:{proteome_id})'
    r = requests.get(url)
    r.raise_for_status()
    with open(f'data/{self.taxon_id}/proteome.fasta', 'w') as f:
      f.write(r.text)

  def _get_proteome_with_most_proteins(self):
    """
    Between UniProt proteome ties, get the proteome with the most proteins.
    """
    # get the index of the row with the most proteins using proteinCount
    proteome_id = self.proteome_list.iloc[self.proteome_list['proteinCount'].idxmax()]['upid']
    proteome_taxon = self.proteome_list.iloc[self.proteome_list['proteinCount'].idxmax()]['taxonomy']

    # write the proteome to a file
    self._get_proteome_to_fasta(proteome_id)

    return proteome_id, proteome_taxon


def main():
  import argparse

  parser = argparse.ArgumentParser()

  parser.add_argument('-t', '--taxon_id', required=True, help='Taxon ID for the species to pull data for.')
  
  args = parser.parse_args()
  taxon_id = args.taxon_id

  Selector = ProteomeSelector(taxon_id)
  print(f'Number of candidate proteomes: {Selector.num_of_proteomes}')
  proteome_data = Selector.select_proteome()
  Selector.proteome_to_csv()

  print(f'Proteome ID: {proteome_data[0]}')
  print(f'Proteome taxon: {int(proteome_data[1])}')
  print(f'Proteome type: {proteome_data[2]}')

if __name__ == '__main__':
  main()