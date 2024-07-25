from fastapi import APIRouter, HTTPException
from services.shared import smiles_to_svg

router = APIRouter()

@router.get("/smiles/draw", tags=['Shared'])
async def draw_smiles(smiles: list[str]):
    return smiles_to_svg(smiles)

# TODO: implement the function
# @router.get("/chemicals/info", tags=['Shared'])
# async def get_chemical_info(ids: list[str]):
    # return None