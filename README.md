# UniProtPy

Python library that interfaces with UniProt API.

For something like [openvax/pyensembl](https://github.com/openvax/pyensembl) with UniProt.

The REST API has changed as of 2022. Many of the ways to extract data from UniProt is now different and there isn't a clean way to interface with it.


### Goals
1. Allow users to pull any kind of data from UniProt.
2. Store and query large data using a local database.
3. Manipulate and output data in many standard formats.

### Installation

```bash
pip install uniprotpy
```

### Getting a proteome for a species
```bash
uniprotpy get-best-proteome --taxon-id 9606
```

### TODO

- Retrieve individual entries in all supported formats.
- Get metadata (protein ID, name, gene, # of isoforms, etc.) for entries.
- Retrieve proteomes via proteome ID or select "best" proteome based on taxon ID.
- Query proteomes for a protein by ID, name, seq, or peptide unit.