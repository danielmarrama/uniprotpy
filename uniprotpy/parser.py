import re
from Bio import SeqIO


def parse_proteome(proteome_file) -> dict:
    """Parse out a proteome FASTA file and return a protein dictionary    
    Args:
        proteome_file: path to a proteome file in FASTA format.
        
    Returns:
        A dictionary mapping protein IDs to keyword-value pairs."""
    proteome_dict = {}
    for record in SeqIO.parse(proteome_file, "fasta"):
        result = parse_protein_string(record.description)
        protein_id = result["protein_id"]
        proteome_dict[protein_id] = parse_protein_string(record.description)
        proteome_dict[protein_id]["sequence"] = str(record.seq)
    return proteome_dict


def parse_protein_string(s):
    pattern = r'^tr\|(\w+)\|(\w+)_\w+\s+(.+) OS=(.+) OX=(\d+) GN=(\w+) PE=(\d+) SV=(\d+)(?: GP=(\d+))?$'
    match = re.match(pattern, s)

    if match:
        return {
            "protein_id": match.group(2),
            "protein_name": match.group(3),
            "species": match.group(4),
            "taxon_id": match.group(5),
            "gene": match.group(6),
            "pe_level": match.group(7),
            "sequence_version": match.group(8) if match.group(8) else "",
            "gene_priority": match.group(9) if match.group(9) else "",
        }
    else:
        raise ValueError(f"Could not parse protein string: {s}")

