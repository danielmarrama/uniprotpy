import argparse

from .helpers import proteome_to_fasta
from .proteome_selector import ProteomeSelector

def parse_arguments():

  parser = argparse.ArgumentParser(description='UniProtPy')

  parser.add_argument(
    '-t', '--taxon_id', help='NCBI taxonomy identifier.'
  )
  parser.add_argument(
    '-p', '--proteome_id', help='UniProt proteome identifier.'
  )
  parser.add_argument(
    '-s', '--store',
    choices=(
      'sql',
      'fasta',
      'csv'
    ),
    default='sql',
    help='\"sql\" will store the proteome in a SQLite database. '
      '\"fasta\" will store the proteome in a FASTA file. '
      '\"csv\" will store the proteome in a CSV file.'
  )

  args = parser.parse_args()

  return args


def run():
  args = parse_arguments()
  
  taxon_id = args.taxon_id
  proteome_id = args.proteome_id
  store = args.store

  if taxon_id:
    assert proteome_id is None, 'Proteome ID cannot be provided when taxon ID is provided.'
    ProteomeSelector(taxon_id).select_best_proteome()
  
  else:
    assert proteome_id is not None, 'Proteome ID must be provided when taxon ID is not provided.'
    assert 'UP' in proteome_id, 'Proteome ID must be a UniProt proteome ID.'
    proteome_to_fasta(proteome_id, compress=False)
