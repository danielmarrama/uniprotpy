from sqlalchemy import Boolean, Column, Float, ForeignKey, Integer, String
from sqlalchemy.orm import relationship, declarative_base

Base = declarative_base()

class UniprotEntry(Base):
    __tablename__ = 'uniprot_entry'
    protein_id = Column(String, primary_key=True)
    protein_name = Column(String)
    species = Column(String)
    taxon_id = Column(Integer)
    gene = Column(String)
    pe_level = Column(Integer)
    sequence_version = Column(Integer)
    gene_priority = Column(Boolean)
    sequence = Column(String)

    def dict(self):
        return {
            "protein_id": self.protein_id,
            "protein_name": self.protein_name,
            "species": self.species,
            "taxon_id": self.taxon_id,
            "gene": self.gene,
            "pe_level": self.pe_level,
            "sequence_version": self.sequence_version,
            "gene_priority": self.gene_priority,
            "sequence": self.sequence
        }