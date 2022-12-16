#!/usr/bin/env python3

import re
import pandas as pd

from Bio import SeqIO

data = []
for protein in list(SeqIO.parse('9606.fasta', 'fasta')):
  sp_or_tr = protein_id = protein.description.split('|')[0]
  protein_id = protein.description.split('|')[1]
  protein_name = re.search(' (.*) OS', protein.description).group(1)
  species = re.search('OS=(.*) OX=', protein.description).group(1)
  taxon = 9606
  try:
    gene = re.search('GN=(.*?) ', protein.description).group(1)
  except AttributeError:
    gene = ''
  try:
    pe_level = int(re.search('PE=(.*?) ', protein.description).group(1))
  except AttributeError:
    pe_level = ''

  seq = protein.seq

  data.append([sp_or_tr, protein_id, protein_name, species, taxon, gene, pe_level, str(seq)])

columns = ['sp_or_tr', 'protein_id', 'protein_name', 'species', 'taxon', 'gene', 'pe_level', 'seq']
df = pd.DataFrame(data, columns=columns)

df.to_csv('9606.csv')