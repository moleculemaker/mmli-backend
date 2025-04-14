import logging
from Bio import UniProt
from typing import TypedDict

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
    structure_url: str
    
class UniprotLookUpByECResult(TypedDict):
    uniprot_id: str
    structure_url: str

UniprotResultDict = dict[str, UniprotRawResultWithStructure | None]

class UniprotService:
    
    structure_url = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb'
    
    @staticmethod
    def get_info_by_accessions(accessions: list[str]) -> UniprotResultDict:
        '''
        Get information about proteins by their accession numbers.
        Uniprot search will return 400 if any of the accession is invalid.

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
            result['structure_url'] = UniprotService.structure_url.format(result['primaryAccession'])
        
        return {result['primaryAccession']: result for result in results}
    
    @staticmethod
    def find_uniprot_ids_for_ec_numbers(ec_number_inputs: list[str], limit: int = 10) -> dict[str, list[UniprotLookUpByECResult]]:
        """
        Get a representative UniProt ID for a given EC number.
        
        Args:
            ec_number (str): The EC number to search for (e.g., "1.1.1.1")
            limit (int): The maximum number of results to return max 50
            
        Returns:
            Optional[str]: The UniProt ID if found, None otherwise
        """
        ec_numbers = list(set(ec_number_inputs))
        limit = min(limit, 50)
        
        ret_val = {}
        for ec_number in ec_numbers:
            if not ec_number in ret_val:
                ret_val[ec_number] = []
                
            raw_results = UniProt.search(
                query=f'ec:{ec_number}',
                fields=['sequence'],
            )[:limit]

            results = list(raw_results)
            
            ret_val[ec_number].extend({
                'uniprot_id': result['primaryAccession'],
                'sequence': result['sequence']['value'],
                'structure_url': UniprotService.structure_url.format(result['primaryAccession'])
            } for result in results)
        
        return ret_val
    
    @staticmethod
    def get_ec_numbers_for_uniprot_ids(uniprot_ids: list[str]) -> dict[str, list[str]]:
        """
        Get EC numbers for a list of UniProt IDs.
        Uniprot search will return 400 if any of the accession is invalid.
        
        Args:
            uniprot_ids (list[str]): List of UniProt IDs

        Returns:
            dict[uniprot_id, list[str]]: Dictionary of UniProt IDs and their EC numbers
        """
        ret_val = {}
        for uniprot_id in uniprot_ids:
            raw_results = UniProt.search(
                query=f'accession:{uniprot_id}',
                fields=['ec'],
            )
            
            results = list(raw_results)
            ec_numbers = set()
            
            # Extract EC numbers from both recommended and alternative names
            for result in results:
                protein_desc = result['proteinDescription']
                
                # Get EC numbers from recommended name
                if 'recommendedName' in protein_desc and 'ecNumbers' in protein_desc['recommendedName']:
                    for ec in protein_desc['recommendedName']['ecNumbers']:
                        if 'value' in ec:
                            ec_numbers.add(ec['value'])
            
            ret_val[uniprot_id] = list(ec_numbers)
        
        return ret_val