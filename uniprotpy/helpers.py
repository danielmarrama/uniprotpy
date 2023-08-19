import requests
import pandas as pd
from pathlib import Path


def get_proteome(proteome_id: str, compress: bool = False) -> None:
    """Get the FASTA file for a proteome from UniProt API.
    
    Args:
      proteome_id: UniProt proteome identifier.
      compress: Whether to download the FASTA file as compressed."""
    compress = 'true' if compress else 'false'
    base_url = 'https://rest.uniprot.org/uniprotkb'
    full_url = f'{base_url}/stream?compressed={compress}&format=fasta&includeIsoform=true&query=(proteome:{proteome_id})'
  
    r = requests.get(full_url)
    r.raise_for_status()
    
    return r.text


def proteome_to_fasta(proteome: str, output_dir: Path) -> None:
    pass

    
def proteome_to_csv(proteome_dict: dict, output_dir: Path) -> None:
    df = pd.DataFrame.from_dict(proteome_dict, orient='index').reset_index(drop=True)
    df.to_csv(output_dir, index=False)


def proteome_to_tsv(proteome_dict: dict, output_dir: Path) -> None:
    df = pd.DataFrame.from_dict(proteome_dict, orient='index').reset_index(drop=True)
    df.to_csv('/home/dan/Desktop/proteome.tsv', sep='\t', index=False)