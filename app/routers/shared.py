from fastapi import APIRouter, HTTPException
from services.shared import draw_chemical_svg

router = APIRouter()

@router.get("/smiles/draw", tags=['Shared'])
async def draw_smiles(smiles: str):
    if type(smiles) != str: 
        raise HTTPException(status_code=400, detail=f"Input must be a single SMILES string. Got type: `{type(smiles)}` with SMILES = `{smiles}`")
    
    try:
        return draw_chemical_svg(smiles)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# TODO: implement the function
# @router.get("/chemicals/info", tags=['Shared'])
# async def get_chemical_info(ids: list[str]):
    # return None