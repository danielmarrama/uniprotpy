#!/usr/bin/env python3

import pandas as pd
import json
import requests
import gzip

directory = '/run/media/dan/HDD/All/Projects/PEPMatch/Proteomes/fasta/test_prod/'

def get_proteomes(taxon_dict):
    '''
    Extracts taxons with proteomes from the UniProt FTP site.
    If get_canonical: get 1 protein per gene proteome as well as total proteome.
    '''

    url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes'

    for taxon_id, proteome_id in taxon_dict.items():
        # check Eukaryota first
        r = requests.get(url + '/Eukaryota/' + proteome_id + '/' + proteome_id + '_' + taxon_id + '.fasta.gz', stream=True)
        
        # then Archaea
        if r.status_code == 404:
            r = requests.get(url + '/Archaea/' + proteome_id + '/' + proteome_id + '_' + taxon_id + '.fasta.gz', stream=True)
            
            # then Bacteria
            if r.status_code == 404:
                r = requests.get(url + '/Bacteria/' + proteome_id + '/' + proteome_id + '_' + taxon_id + '.fasta.gz', stream=True)
                
                # then Viruses
                if r.status_code == 404:
                    r = requests.get(url + '/Viruses/' + proteome_id + '/' + proteome_id + '_' + taxon_id + '.fasta.gz', stream=True)
                    
                    # if found nowhere, move on
                    if r.status_code == 404:
                        print('Did not retrieve:', taxon_id, 'Proteome ID:', proteome_id)
                        continue
        
        with open(directory + taxon_id + '.fasta', 'wb') as f1:
            f1.write(gzip.open(r.raw, 'rb').read())

    return 0


def get_protein_entries(taxons):
    '''
    Get collection of proteins for taxons with no proteome in UniProt.
    '''

    for taxon in taxons:
        r = requests.get('https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A' + taxon + '%29')
        with open(directory + taxon + '.fasta', 'w') as f:
            f.write(r.text)

if __name__ == '__main__':

    # get taxons needed for the PEPMatch server
    with open('taxon_id_data.json', 'r') as f:
        j = json.load(f)
    taxons = pd.Series(j.keys())

    # read in the table from UniProt (https://www.uniprot.org/proteomes?query=*)
    df = pd.read_csv('proteome_ids.tsv', sep='\t')

    # get taxons with proteomes (those that have IDs)
    taxons_with_proteomes = list(taxons[taxons.isin(df['Organism Id'].astype(str))])

    # get proteome ID with most proteins
    df.drop_duplicates(subset=['Organism Id'], inplace=True)
    # idx = df.groupby(['Organism Id'])['Protein count'].transform(max) == df['Protein count']
    # df_temp = df[idx].drop_duplicates(subset=['Organism Id'])
    taxons_with_proteomes_dict = dict(zip(df['Organism Id'].astype(str), df['Proteome Id'].astype(str)))
    
    # now get those proteomes from UniProt FTP server
    get_proteomes(taxons_with_proteomes_dict)
    

    # get collection of proteins for taxons that have no proteome
    taxons_without_proteomes = list(taxons[~taxons.isin(df['Organism Id'].astype(str))])
    get_protein_entries(taxons_without_proteomes)