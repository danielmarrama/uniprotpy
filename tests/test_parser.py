import pytest
from unittest.mock import MagicMock, patch
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

from uniprotpy.parser import parse_proteome


test_result_expected = {
    "protein_id": "A0A075B6G3",
    "protein_name": "Dystrophin",
    "species": "Homo sapiens",
    "taxon_id": "9606",
    "gene": "DMD",
    "pe_level": "1",
    "sequence_version": "1",
    "gene_priority": "0",
    "sequence": "MLWWEEVEDCYEREDVQKKTFTKWVNAQFSKFGKQHIENLFSDLQDGRRLLDLLEGLTGQ",
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
        assert result[test_result_expected["protein_id"]] == test_result_expected
