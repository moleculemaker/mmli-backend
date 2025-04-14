import logging
from Bio import UniProt
from typing import TypedDict, Any

log = logging.getLogger(__name__)

class UniprotException(BaseException):
    type: str
    message: str
    
    def __init__(self, type: str, message: str):
        self.type = type
        self.message = message
        

class UniprotOrgansimEntry(TypedDict):
    scientificName: str
    commonName: str
    taxonId: int
    lineage: list[str]
    
class UniprotNameEntry(TypedDict):
    value: str
    
class UniprotNamesEntry(TypedDict):
    fullName: UniprotNameEntry
    shortNames: list[UniprotNameEntry]
    
class UniprotProteinDescription(TypedDict):
    recommendedName: UniprotNamesEntry
    alternativeNames: list[UniprotNamesEntry]
    flag: str

class UniprotGeneEntry(TypedDict):
    geneName: UniprotNameEntry
    
class UniprotExtraAttributesEntry(TypedDict):
    uniParcId: str
    
class UniprotSequenceEntry(TypedDict):
    value: str
    length: int
    molWeight: int
    crc64: str
    md5: str
    
class UniprotRawResult(TypedDict):
    entryType: str
    primaryAccession: str
    organism: UniprotOrgansimEntry
    proteinDescription: UniprotProteinDescription
    genes: list[UniprotGeneEntry]
    sequence: UniprotSequenceEntry
    extraAttributes: UniprotExtraAttributesEntry
    
class UniprotRawResultWithStructure(UniprotRawResult):
    url: str

UniprotResultDict = dict[str, UniprotRawResultWithStructure | None]

class UniprotService:
    
    @staticmethod
    def get_info_by_accessions(accessions: list[str]) -> UniprotResultDict:
        '''
        Get information about proteins by their accession numbers.
        Uniprot search will return 400 if any of the accessions are invalid.

        Parameters:
            accessions (list[str]): List of UniProt accession numbers

        Returns:
            dict[accession, UniprotRawResult | None]: Dictionary of protein information dictionaries
        '''
        rawResults = UniProt.search(
            query=' OR '.join(f'accession:{accession}' for accession in accessions), 
            fields=['gene_names', 'organism_name', 'lineage', 'protein_name', 'sequence'], 
            batch_size=10
        )
        
        results = list(rawResults)
        for result in results:
            result['url'] = f'https://alphafold.ebi.ac.uk/files/AF-{result["primaryAccession"]}-F1-model_v4.pdb'
        
        return {result['primaryAccession']: result for result in results}