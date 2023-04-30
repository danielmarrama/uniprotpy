import os
import tempfile
import pytest

from uniprotpy.database import UniprotDatabase

@pytest.fixture(scope='module')
def temp_db_url():
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.close()
        yield f"sqlite:///{temp_file.name}"
        os.unlink(temp_file.name)

@pytest.fixture(scope='module')
def database(temp_db_url):
    db = UniprotDatabase(database_path=temp_db_url)
    db._init_sqlite()
    return db


test_data = {
    "protein_id": "|A0A075B6G3|",
    "protein_name": " Dystrophin OS",
    "species": "OS=Homo sapiens OX",
    "taxon_id": "OX=9606 ",
    "gene": "GN=DMD ",
    "pe_level": "PE=1 ",
    "sequence_version": "",
    "gene_priority": "",
    "sequence": "MLWWEEVEDCYEREDVQKKTFTKWVNAQFSKFGKQHIENLFSDLQDGRRLLDLLEGLTGQ",
}

def test_add_protein_entry(database):
    database.add(test_data)

def test_get_protein_entry(database):
    result = database.get(test_data["protein_id"]).dict()
    assert result["protein_id"] == "|A0A075B6G3|"
