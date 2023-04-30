from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker
from uniprotpy.models import UniprotEntry, Base

class UniprotDatabase():
    def __init__(self, species=None, proteome_id=None, database_path=None):
        self.species = species
        self.proteome_id = proteome_id
        self.database_path = database_path
        self.engine = create_engine(self.database_path)
        self.session = sessionmaker(bind=self.engine)
        self._init_sqlite()

    def _init_sqlite(self):
        """Initialize a sqlite database."""
        if not inspect(self.engine).has_table("uniprot_entry"):
            Base.metadata.create_all(self.engine)

    def add(self, protein):
        """Given a dictionary containing a uniprot entry, add it to the database.

        Args:
            protein (dict): Dictionary containing a uniprot entry.
        """
        session = self.session()
        session.add(UniprotEntry(**protein))
        session.commit()
        session.close()

    def get(self, protein_id):
        """Given a protein ID, return the corresponding uniprot entry.

        Args:
            protein_id (str): Protein ID.
        """
        session = self.session()
        result = session.query(UniprotEntry).filter_by(protein_id=protein_id).first()
        session.close()
        return result
    
    
    def list(self):
        """Return a list of all protein IDs in the database."""
        session = self.session()
        result = session.query(UniprotEntry).all()
        session.close()
        return result