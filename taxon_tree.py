#!/usr/bin/env python3

import requests
import pandas as pd
import io

BASE_URL = 'https://rest.uniprot.org/taxonomy/stream'


def get_taxon_children(taxon_id):
  url = f'{BASE_URL}?fields=id%2Cscientific_name&format=tsv&query=%28parent%3A{taxon_id}%29'
  response = requests.get(url)

  try:
    taxon_children_df = pd.read_csv(io.StringIO(response.text), sep='\t')
    children = taxon_children_df.to_dict(orient='records')
  except ValueError:
    children = []

  return children


def traverse_tree(taxon_id):
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


def main():
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--taxon_id', type=int, help='Taxon ID')
  parser.add_argument('-p', '--pickle', action='store_true', help='Pickle output.')

  args = parser.parse_args()
  taxon_id = args.taxon_id

  taxonomy_tree = traverse_tree(taxon_id)

  if args.pickle:
    import pickle
    with open('taxonomy_tree.pkl', 'wb') as f:
      pickle.dump(taxonomy_tree, f)
  
  print('Taxon tree:')
  print(taxonomy_tree)

if __name__ == '__main__':
  main()