from fastapi import APIRouter, HTTPException
from services.shared import draw_chemical_svg

router = APIRouter()

@router.get("/smiles/draw", tags=['Shared'])
async def draw_smiles(smiles: list[str]):
    return draw_chemical_svg(smiles)

# TODO: implement the function
# @router.get("/chemicals/info", tags=['Shared'])
# async def get_chemical_info(ids: list[str]):
    # return None