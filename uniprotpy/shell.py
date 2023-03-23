import argparse

from .select_best_proteome import ProteomeSelector

def parse_arguments():
  parser = argparse.ArgumentParser(description='UniProtPy')
  parser.add_argument('-t', '--taxon_id', required=True)
  parser.add_argument(
    'action',
    choices=(
      'download',
    ),
    help='\"download\" will download the best proteome from UniProt for a given '
         'taxon ID of a species.'
  )
  args = parser.parse_args()

  return args

def run():
  args = parse_arguments()
  taxon_id = args.taxon_id
  if args.action == 'download':
    ProteomeSelector(taxon_id).select_proteome()
  else:
    raise ValueError('Invalid action.')