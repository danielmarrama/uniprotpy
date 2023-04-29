import pytest
from unittest.mock import MagicMock, patch
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

from uniprotpy.parser import parse_proteome

test_result_expected = {
    "|A0A075B6G3|": {
        "protein_name": " Dystrophin OS",
        "species": "OS=Homo sapiens OX",
        "taxon_id": "OX=9606 ",
        "gene": "GN=DMD ",
        "pe_level": "PE=1 ",
        "sequence_version": "",
        "gene_priority": "",
        "sequence": "MLWWEEVEDCYEREDVQKKTFTKWVNAQFSKFGKQHIENLFSDLQDGRRLLDLLEGLTGQ",
    }
}


def test_parse_proteome():
    mock_seq_record = SeqRecord(
        Seq("MLWWEEVEDCYEREDVQKKTFTKWVNAQFSKFGKQHIENLFSDLQDGRRLLDLLEGLTGQ"),
        id="tr|A0A075B6G3|A0A075B6G3_HUMAN",
        description="tr|A0A075B6G3|A0A075B6G3_HUMAN Dystrophin OS=Homo sapiens OX=9606 GN=DMD PE=1 SV=1",
    )
    mock_seqio_parse = MagicMock(return_value=[mock_seq_record])

    with patch("Bio.SeqIO.parse", mock_seqio_parse):
        # Call the parse_proteome function with the mocked objects
        result = parse_proteome("dummy_path")
        assert result == test_result_expected
