import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_proteome(proteome_file) -> dict:
    """Parse out a proteome FASTA file and return a protein dictionary    
    Args:
        proteome_file: path to a proteome file in FASTA format.
        
    Returns:
        A dictionary mapping protein IDs to keyword-value pairs."""
    proteome_dict = {}
    for record in SeqIO.parse(proteome_file, "fasta"):
        protein_data = parse_protein_record(record)
        protein_id = protein_data["protein_id"]
        proteome_dict[protein_id] = protein_data
        proteome_dict[protein_id]["sequence"] = str(record.seq)
    return proteome_dict


def parse_protein_record(record: SeqRecord) -> dict:
    """Parse a record from a FASTA file and return a list of the protein data.
    Args:
        record: Biopython SeqRecord of the protein entry.
        
    Returns:
        A list of the protein data in the order of the regexes below."""
    regexes = {
        'protein_id': re.compile(r"\|([^|]*)\|"),     # between | and |
        'protein_name': re.compile(r"\s(.+?)\sOS"),   # between space and space before OS
        'species': re.compile(r"OS=(.+?)\sOX"),       # between OS= and space before OX
        'taxon_id': re.compile(r"OX=(.+?)(\s|$)"),         # between OX= and space
        'gene': re.compile(r"GN=(.+?)(\s|$)"),             # between GN= and space
        'pe_level': re.compile(r"PE=(.+?)(\s|$)"),         # between PE= and space
        'sequence_version': re.compile(r"SV=(.+?)(\s|$)"), # between SV= and space
        'gene_priority': re.compile(r"GP=(.+?)(\s|$)"),    # between GP= and space
    }
    protein_data = {}
    for key in regexes: # loop through compiled regexes to extract protein data
        match = regexes[key].search(str(record.description))
        
        if match:
            protein_data[key] = match.group(1)
        else:
            if key == 'protein_id':
                protein_data[key] = str(record.id) # get record.id from FASTA header instead
            elif key == 'sequence_version':
                protein_data[key] = '1'
            elif key in ['pe_level', 'gene_priority']:
                protein_data[key] = '0' # zeros for integer columns
            else:
                protein_data[key] = ''  # empty strings for string columns

    return protein_data
