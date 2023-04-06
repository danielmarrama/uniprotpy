#!/usr/bin/env python3

import requests


def proteome_to_fasta(proteome_id, compress=False):
  """
  Get the FASTA file for a proteome from UniProt API.
  """
  compress = 'true' if compress else 'false'
  url = f'https://rest.uniprot.org/uniprotkb/stream?compressed={compress}&format=fasta&includeIsoform=true&query=(proteome:{proteome_id})'
  r = requests.get(url)
  r.raise_for_status()
  with open(f'{proteome_id}.fasta', 'w') as f:
    f.write(r.text)