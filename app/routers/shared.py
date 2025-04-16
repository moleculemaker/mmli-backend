import traceback
from typing import Optional, TypeVar, TypedDict
from fastapi import APIRouter, Body, HTTPException

from services.shared import (
    ChemInfoAlgorithms,
    ChemInfoAlgorithmsConfig,
    convert_to_rdkit_mols, 
    draw_chemical_svg, 
    rdkit_get_cheminfo_algorithms_results, 
)
from config import get_logger

router = APIRouter()

log = get_logger(__name__)

T = TypeVar('T')
DictResponse = dict[str, T]

class FragmentMatchResponseData(TypedDict):
    matches: list[list[int]]
    # svg: str
    
class ChemInfoResponse(TypedDict):
    tanimoto: Optional[dict[str, float]]
    fragment: Optional[dict[str, FragmentMatchResponseData]]
    mcs: Optional[dict[str, float]]
    errors: Optional[dict[str, str | dict[str, str]]]

@router.get("/smiles/draw", tags=['Shared'])
async def draw_smiles(smiles: str):
    if type(smiles) != str: 
        raise HTTPException(status_code=400, detail=f"Input must be a single SMILES string. Got type: `{type(smiles)}` with SMILES = `{smiles}`")
    
    try:
        return draw_chemical_svg(smiles)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    

@router.post("/rdkit/cheminfo-algorithms", 
            tags=['Shared'], 
            response_model=ChemInfoResponse)
async def get_cheminfo_algorithms_results(
    algorithms: list[ChemInfoAlgorithms] = Body(...),
    query_smiles: str = Body("COC1=C(C=CC(=C1)CCN)O"),
    smiles_list: list[str] = Body(["COc1ccc2c(c1OC)C(O)O[C@@H]2[C@H]1c2c(cc3c(c2OC)OCO3)CCN1C"]),
    algorithm_config: Optional[ChemInfoAlgorithmsConfig] = Body(ChemInfoAlgorithmsConfig()),
):
    """
    Run cheminformatics algorithms on a list of SMILES strings against a query SMILES string.
    """
    # Convert SMILES to molecules
    smi_mol_dict, mol_smi_dict = convert_to_rdkit_mols(smiles_list + [query_smiles])
    
    # Validate query SMILES
    query_mol = smi_mol_dict[query_smiles]
    if query_mol is None:
        raise HTTPException(
            status_code=400, 
            detail=f"fail to process query SMILES string: {query_smiles}"
        )
        
    # Format results for response
    response = {
        "errors": {},
        "tanimoto": None,
        "fragment": None,
        "mcs": None,
    }
    
    for smiles, mol in smi_mol_dict.items():
        if mol is None:
            response["errors"][smiles] = f"fail to process SMILES string: {smiles}"
    
    try:
        # Get results from cheminformatics algorithms
        results = rdkit_get_cheminfo_algorithms_results(
            algorithms, 
            query_mol, 
            [ v for k, v in smi_mol_dict.items() if k != query_smiles and v is not None ],
            algorithm_config
        )
        
        if "tanimoto" in algorithms:
            response["tanimoto"] = { 
                mol_smi_dict[result["mol"]]: result["score"] for result in results["tanimoto"] if "error" not in result
            }
            response["errors"]["tanimoto"] = {
                mol_smi_dict[result["mol"]]: result["error"] for result in results["tanimoto"] if "error" in result
            }
            
        if "fragment" in algorithms:
            response["fragment"] = { 
                mol_smi_dict[result["mol"]]: {
                    "matches": result["matches"], 
                    # "svg": result["svg"] 
                } for result in results["fragment"] if "error" not in result
            }
            response["errors"]["fragment"] = {
                mol_smi_dict[result["mol"]]: result["error"] for result in results["fragment"] if "error" in result
            }
            
        if "mcs" in algorithms:
            response["mcs"] = { 
                mol_smi_dict[result["mol"]]: result["score"] for result in results["mcs"] if "error" not in result
            }
            response["errors"]["mcs"] = {
                mol_smi_dict[result["mol"]]: result["error"] for result in results["mcs"] if "error" in result
            }
            
        return response
        
    except Exception as e:
        log.error(f"Error processing cheminformatics algorithms: {str(e)}, {traceback.format_exc()}")
        raise HTTPException(
            status_code=500,
            detail=f"Error processing cheminformatics algorithms: {str(e)}"
        )
  
# it's available but currently not used  
# @router.post(
#     "/uniprot/get-info-by-accessions", 
#     tags=['Shared'],
#     response_model=UniprotResultDict,
#     description="""
#     Get information about proteins by their accession numbers.
#     Uniprot search will return 400 if any of the accessions are invalid.
#     Returns a dictionary with the accession numbers as keys and the protein information as values.
#     """,
# )
# def get_info_by_accessions(accessions: list[str]):
#     try:
#         return UniprotService.get_info_by_accessions(accessions)
#     except Exception as e:
#         raise HTTPException(status_code=500, detail=str(e))
    
# @router.post(
#     "/kegg/get-info-by-ec-numbers",
#     tags=['Shared'],
#     response_model=KeggResultDict,
#     description="""
#     Get information about enzymes by their EC numbers.
#     """,
# )
# def get_info_by_ec_numbers(ec_numbers: list[str]):
#     try:
#         return KeggService.get_info_by_ec_numbers(ec_numbers)
#     except Exception as e:
#         raise HTTPException(status_code=500, detail=str(e))
    
# TODO: implement the function
# @router.get("/chemicals/info", tags=['Shared'])
# async def get_chemical_info(ids: list[str]):
    # return None