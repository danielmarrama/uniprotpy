import re
from Bio import SeqIO

regexes = {
        'protein_id': re.compile(r"\|([^|]*)\|"),      # between | and |
        'protein_name': re.compile(r"\s(.+?)\sOS"),    # between first space and space before OS
        'species': re.compile(r"OS=(.+?)\sOX"),        # between OS= and space before OX (species can have spaces)
        'taxon_id': re.compile(r"OX=(.+?)\s"),         # between OX= and space
        'gene': re.compile(r"GN=(.+?)\s"),             # between GN= and space
        'pe_level': re.compile(r"PE=(.+?)\s"),         # between PE= and space
        'sequence_version': re.compile(r"SV=(.+?)\s"), # between SV= and space
        'gene_priority': re.compile(r"GP=(.+?)\s"),    # between GP= and space
}

def find_pattern(pattern, string):
    match = re.search(pattern, string)
    if match:
        return match.group()
    else:
        return ""

def parse_proteome(proteome_file):
    proteome_dict = {}
    for record in SeqIO.parse(proteome_file, "fasta"):
        protein_id = find_pattern(regexes['protein_id'], record.description)
        protein_name = find_pattern(regexes['protein_name'], record.description)
        species = find_pattern(regexes['species'], record.description)
        taxon_id = find_pattern(regexes['taxon_id'], record.description)
        gene = find_pattern(regexes['gene'], record.description)
        pe_level = find_pattern(regexes['pe_level'], record.description)
        sequence_version = find_pattern(regexes['sequence_version'], record.description)
        gene_priority = find_pattern(regexes['gene_priority'], record.description)
        
        proteome_dict[protein_id] = {
            "protein_name": protein_name,
            "species": species,
            "taxon_id": taxon_id,
            "gene": gene,
            "pe_level": pe_level,
            "sequence_version": sequence_version,
            "gene_priority": gene_priority,
            "sequence": str(record.seq)
        }
    return proteome_dict