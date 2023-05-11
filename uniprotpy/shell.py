import argparse

from .helpers import proteome_to_fasta
from .proteome_selector import ProteomeSelector

def parse_arguments():
  parser = argparse.ArgumentParser(description='UniProtPy')
  parser.add_argument('-t', '--taxon_id', help='NCBI taxonomy identifier.')
  parser.add_argument('-p', '--proteome_id', help='UniProt proteome identifier.')
  parser.add_argument(
    'action',
    choices=(
      'get-proteome',
      'get-best-proteome',
    ),
    help='\"get-proteome\" will download a proteome from UniProt using a proteome ID. ' 
         '\"get-best-proteome\" will download the best proteome from UniProt for a given '
         'taxon ID of species rank or below.'
  )
  args = parser.parse_args()

  return args

def run():
  args = parse_arguments()
  taxon_id = args.taxon_id
  proteome_id = args.proteome_id

  if args.action == 'get-best-proteome':
    assert taxon_id is not None
    ProteomeSelector(taxon_id).select_proteome()
  elif args.action == 'get-proteome':
    assert proteome_id is not None
    proteome_to_fasta(proteome_id, compress=False)
  else:
    raise ValueError('Invalid action.')