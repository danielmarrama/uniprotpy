#!/usr/bin/env python3

import requests
import pandas as pd
import io

BASE_URL = 'https://rest.uniprot.org/taxonomy/stream'


def get_taxon_children(taxon_id):
  """Retrives all children taxa for a given taxon ID."""
  url = f'{BASE_URL}?fields=id%2Cscientific_name&format=tsv&query=%28parent%3A{taxon_id}%29'
  response = requests.get(url)

  try:
    taxon_children_df = pd.read_csv(io.StringIO(response.text), sep='\t')
    children = taxon_children_df.to_dict(orient='records')
  except ValueError:
    children = []

  return children

def traverse_tree(taxon_id):
  """
  Recursively traverse the taxonomy tree for a given taxon ID.

  It goes down to each leaf node and returns a dict of the form:
  {taxon_id: [child1, {child2: subchild1}, ...]}
  
  Taxon IDs that have no children are standalone integers; anything
  else is a list of dicts.
  """
  children = get_taxon_children(taxon_id)

  if not children:
    return taxon_id

  child_results = []
  for child in children:
    child_id = child['Taxon Id']
    child_result = traverse_tree(child_id)

    if isinstance(child_result, int):
      child_results.append(child_result)
    else:
      child_results.append(child_result)

  return {taxon_id: child_results}

def create_taxon_tree(tree, prefix="", file=None):
  """
  Create a taxonomy tree structure from the output of traverse_tree. 

  It looks like this for the homo genus (9605):
  └──9605
      ├──9606
      │   ├──63221
      │   └──741158
      ├──1425170
      ├──2665952
      │   └──2665953
      └──2813598
          └──2813599
  """
  if isinstance(tree, list):
    for i, item in enumerate(tree):
      is_last = i == len(tree) - 1
      if isinstance(item, dict):
        for k, v in item.items():
          print(f"{prefix}{'└──' if is_last else '├──'}{k}", file=file)
          create_taxon_tree(v, f"{prefix}{'    ' if is_last else '│   '}", file)
      else:
        print(f"{prefix}{'└──' if is_last else '├──'}{item}", file=file)
  
  elif isinstance(tree, dict):
    for i, (key, value) in enumerate(tree.items()):
      is_last = i == len(tree) - 1
      print(f"{prefix}{'└──' if is_last else '├──'}{key}", file=file)
      create_taxon_tree(value, f"{prefix}{'    ' if is_last else '│   '}", file)

def main():
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--taxon_id', required=True, type=int, help='Taxon ID')
  parser.add_argument('-p', '--pickle', action='store_true', help='Pickle tree output.')
  parser.add_argument('-o', '--output', action='store_true', help='Output tree file.')

  args = parser.parse_args()
  taxon_id = args.taxon_id

  taxonomy_tree = traverse_tree(taxon_id)

  if args.pickle:
    import pickle
    with open('taxonomy_tree.pkl', 'wb') as f:
      pickle.dump(taxonomy_tree, f)
  
  if args.output:
    with open(f'{taxon_id}_taxonomy_tree.txt', 'w') as f:
      create_taxon_tree(taxonomy_tree, file=f)

  print(f'Taxon tree for {taxon_id}:')
  create_taxon_tree(taxonomy_tree)

if __name__ == '__main__':
  main()