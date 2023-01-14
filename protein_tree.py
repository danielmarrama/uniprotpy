#!/usr/bin/env python3

import re
import pandas as pd

from Bio import SeqIO
from anytree import Node, RenderTree


def add_nodes(nodes, parent, child):
    if parent not in nodes:
        nodes[parent] = Node(parent)  
    if child not in nodes:
        nodes[child] = Node(child)
    nodes[child].parent = nodes[parent]

def create_protein_tree(proteome):
  proteins = list(SeqIO.parse(proteome, 'fasta'))

  data = []
  for protein in proteins:
    try:
      gene = re.search('GN=(.*?) ', protein.description).group(1)
    except AttributeError:
      try:
        gene = re.search('GN=(.*?)$', protein.description).group(1)
      except AttributeError:
        gene = ''
    
    data.append([protein.id.split('|')[0], gene, protein.id.split('|')[1], str(protein.seq)])
  
  df = pd.DataFrame(data, columns=['db', 'gene', 'id', 'seq'])

  nodes = {}
  for parent, child in zip(df['gene'],df['id']):
      add_nodes(nodes, parent, child)

  with open('protein_tree.txt', 'w') as f:
    roots = list(df[~df['gene'].isin(df['id'])]['gene'].unique())
    for root in roots:    
        for pre, _, node in RenderTree(nodes[root]):
            f.write("%s%s" % (pre, node.name))
            f.write('\n')

if __name__ == '__main__':
  proteome = '9606.fasta'
  create_protein_tree(proteome)