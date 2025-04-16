import logging
import time
from typing import Tuple, TypedDict

from Bio.KEGG import REST
from Bio.KEGG import Enzyme

from services.uniprot_service import UniprotLookUpByECResult, UniprotService
from utils.retry import retry

log = logging.getLogger(__name__)

class KeggException(BaseException):
    type: str
    message: str
    
    def __init__(self, type: str, message: str):
        self.type = type
        self.message = message
        
class KeggPathway(TypedDict):
    database: str
    id: str
    description: str
    url: str
    
class KeggGeneEntry(TypedDict):
    id: str
    url: str
    
class KeggGenesByOrganism(TypedDict):
    organism: str
    entries: list[KeggGeneEntry]
    url: str

class KeggRawResult(TypedDict):
    entry: str
    name: list[str]
    classname: list[str]
    sysname: list[str]
    reaction: list[str]
    substrate: list[str]
    product: list[str]
    inhibitor: list[str]
    cofactor: list[str]
    effector: list[str]
    comment: list[str]
    
    pathway: list[KeggPathway]
    
    # Temporary disable genes parsing because there are too many genes
    # genes: list[KeggGenesByOrganism]
    
    # TODO: parse into urls
    disease: list[Tuple[str, str, str]]     # (database, id, disease)
    structures: list[Tuple[str, list[str]]] # (database, list of struct ids)
    dblinks: list[Tuple[str, list[str]]]    # (database, list of db ids)
    
    # get representative uniprot id for each ec number
    uniprots: list[UniprotLookUpByECResult]

KeggResultDict = dict[str, KeggRawResult | None]

class KeggService:
    
    @staticmethod
    @retry(max_retries=3, initial_delay=1.0, max_delay=10.0, backoff_factor=2.0)
    def get_info_by_ec_numbers(ec_numbers_input: list[str]) -> KeggResultDict:
        '''
        Get information about ec numbers. 
        This function is rate limited to 3 requests per second.
        Return 404 if the ec number is not found.
        
        Parameters:
            ec_numbers (list[str]): List of EC numbers

        Returns:
            dict[ec_number, KeggRawResult | None]: Dictionary of EC number information dictionaries
        '''
        BATCH_SIZE = 10 # Bio.KEGG max entries per request
        RATE_LIMIT = 3 # requests per second
        all_records: list[Enzyme.Record] = []
        ec_numbers = list(set(ec_numbers_input))
        uniprot_limit = 3
        
        # Get uniprot ids for ec numbers
        uniprot_ids = UniprotService.find_uniprot_ids_for_ec_numbers(ec_numbers, uniprot_limit)
        
        for i in range(0, len(ec_numbers), BATCH_SIZE):
            batch = ec_numbers[i:i + BATCH_SIZE]
            result_handle = REST.kegg_get(["ec:" + ec_number for ec_number in batch])
            
            records_str = result_handle.read()
            lines = records_str.splitlines()
            parsed_records = list(Enzyme.parse(lines))
            all_records.extend(parsed_records)
            
            # Sleep for 1 second to avoid rate limiting
            if i % RATE_LIMIT == 0 and i != 0:
                time.sleep(1)
            
        result_dict: KeggResultDict = {}
        
        for record in all_records:
            if not record.entry:
                continue
                
            # Extract pathways
            pathways = []
            if hasattr(record, 'pathway'):
                for pathway in record.pathway:
                    if len(pathway) >= 2:
                        pathways.append({
                            'database': pathway[0],
                            'id': pathway[1],
                            'description': pathway[2] if len(pathway) > 2 else '',
                            'url': f"https://www.kegg.jp/pathway/{pathway[1]}"
                        })
            
            # Extract genes
            genes = []
            if hasattr(record, 'genes'):
                for gene_entry in record.genes:
                    if len(gene_entry) >= 2:
                        organism = gene_entry[0]
                        gene_entries = []
                        
                        for gene in gene_entry[1]:
                            gene_entries.append({
                                'id': gene,
                                'url': f"https://www.kegg.jp/entry/{organism}+{gene}"
                            })
                        genes.append({
                            'organism': organism,
                            'entries': gene_entries,
                            'url': f"https://www.kegg.jp/entry/{organism}+{'+'.join(gene_entry[1])}"
                        })
            
            # Create KeggRawResult object
            kegg_result = KeggRawResult(
                entry=record.entry,
                name=record.name if hasattr(record, 'name') else [],
                classname=record.classname if hasattr(record, 'classname') else [],
                sysname=record.sysname if hasattr(record, 'sysname') else [],
                reaction=record.reaction if hasattr(record, 'reaction') else [],
                substrate=record.substrate if hasattr(record, 'substrate') else [],
                product=record.product if hasattr(record, 'product') else [],
                inhibitor=record.inhibitor if hasattr(record, 'inhibitor') else [],
                cofactor=record.cofactor if hasattr(record, 'cofactor') else [],
                effector=record.effector if hasattr(record, 'effector') else [],
                comment=record.comment if hasattr(record, 'comment') else [],
                pathway=pathways,
                # genes=genes,
                disease=record.disease,  # TODO: Implement disease parsing
                structures=record.structures,  # TODO: Implement structures parsing
                dblinks=record.dblinks,  # TODO: Implement dblinks parsing
                uniprots=uniprot_ids[record.entry]
            )
            
            result_dict[record.entry] = kegg_result
            
        return result_dict
    

# Testing
# if __name__ == "__main__":
#     import json
#     result = KeggService.get_info_by_ec_numbers(["1.1.1.1", "1.1.1.2"])
#     with open("kegg_result.json", "w") as f:
#         json.dump(result, f, indent=4)