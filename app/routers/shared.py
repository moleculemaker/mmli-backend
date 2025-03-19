from typing import Dict
from fastapi import APIRouter, Body, HTTPException
from services.shared import draw_chemical_svg, get_similarity_score

router = APIRouter()

@router.get("/smiles/draw", tags=['Shared'])
async def draw_smiles(smiles: str):
    if type(smiles) != str: 
        raise HTTPException(status_code=400, detail=f"Input must be a single SMILES string. Got type: `{type(smiles)}` with SMILES = `{smiles}`")
    
    try:
        return draw_chemical_svg(smiles)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/smiles/similarity", tags=['Shared'], response_model=Dict[str, float | str])
async def get_similarity(
    query_smiles: str = Body(...),
    smiles_list: list[str] = Body(...)
):
    
    try:
        err_code, ret_val = get_similarity_score(query_smiles, smiles_list)
        if err_code != 0:
            raise HTTPException(status_code=400, detail=ret_val)
        else:
            return ret_val

    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
    
# TODO: implement the function
# @router.get("/chemicals/info", tags=['Shared'])
# async def get_chemical_info(ids: list[str]):
    # return None