from fastapi import APIRouter, HTTPException
from services.kegg_service import KeggResultDict, KeggService
from services.shared import draw_chemical_svg
from services.uniprot_service import UniprotService, UniprotResultDict

router = APIRouter()

@router.get("/smiles/draw", tags=['Shared'])
async def draw_smiles(smiles: str):
    if type(smiles) != str: 
        raise HTTPException(status_code=400, detail=f"Input must be a single SMILES string. Got type: `{type(smiles)}` with SMILES = `{smiles}`")
    
    try:
        return draw_chemical_svg(smiles)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
@router.post(
    "/uniprot/get-info-by-accessions", 
    tags=['Shared'],
    response_model=UniprotResultDict,
    description="""
    Get information about proteins by their accession numbers.
    Uniprot search will return 400 if any of the accessions are invalid.
    Returns a dictionary with the accession numbers as keys and the protein information as values.
    """,
)
def get_info_by_accessions(accessions: list[str]):
    try:
        return UniprotService.get_info_by_accessions(accessions)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
@router.post(
    "/kegg/get-info-by-ec-numbers",
    tags=['Shared'],
    response_model=KeggResultDict,
    description="""
    Get information about enzymes by their EC numbers.
    """,
)
def get_info_by_ec_numbers(ec_numbers: list[str]):
    try:
        return KeggService.get_info_by_ec_numbers(ec_numbers)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# TODO: implement the function
# @router.get("/chemicals/info", tags=['Shared'])
# async def get_chemical_info(ids: list[str]):
    # return None