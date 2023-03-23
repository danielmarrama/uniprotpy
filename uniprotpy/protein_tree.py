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

  # get UniProt IDs for one protein per gene proteome
  gp_ids = [str(x.id.split('|')[1]) for x in list(SeqIO.parse(gp_proteome, 'fasta'))]

  data = []
  for protein in proteins:

    # get gene symbol from FASTA file
    try:
      gene = re.search('GN=(.*?) ', protein.description).group(1)
    except AttributeError:

      # sometimes the gene symbol is at the end of the FASTA description
      try:
        gene = re.search('GN=(.*?)$', protein.description).group(1)
      except AttributeError:
        gene = ''

    data.append([protein.id.split('|')[0], gene, protein.id.split('|')[1], str(protein.seq)])

  # put protein tree data into dataframe
  df = pd.DataFrame(data, columns=['db', 'gene', 'id', 'seq'])

  # start tree with nodes - genes as root and UniProt IDs as children
  nodes = {}
  for parent, child in zip(df['gene'],df['id']):
    add_nodes(nodes, parent, child)

  # write the tree into a text file
  with open('protein_tree.txt', 'w') as f:
    roots = list(df[~df['gene'].isin(df['id'])]['gene'].unique())
    for root in roots:         # you can skip this for roots[0], if there is no forest and just 1 tree
      for pre, _, node in RenderTree(nodes[root]):
        if node.name in gp_ids:
          f.write("%s%s*" % (pre, node.name))
          f.write('\n')
        else:
          f.write("%s%s" % (pre, node.name))
          f.write('\n')
