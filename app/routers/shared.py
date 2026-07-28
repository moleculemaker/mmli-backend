from typing import Annotated
from fastapi import APIRouter, File, HTTPException, Query, UploadFile
from services.shared import draw_chemical_svg, is_valid_pdb_file
from services.rdkit_service import RDKitService

router = APIRouter()

@router.get("/smiles/draw", tags=['Shared'])
async def draw_smiles(
    smiles: str,
    highlightAtoms: Annotated[list[int] | None, Query()] = Query(default=None)
):
    if type(smiles) != str:
        raise HTTPException(status_code=400, detail=f"Input must be a single SMILES string. Got type: `{type(smiles)}` with SMILES = `{smiles}`")

    try:
        return draw_chemical_svg(smiles, highlightAtoms=highlightAtoms)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# TODO: implement the function
# @router.get("/chemicals/info", tags=['Shared'])
# async def get_chemical_info(ids: list[str]):
    # return None

@router.get("/smiles/canonicalize", tags=['Shared'])
async def canonicalize_smiles(
    smiles: str
):
    if not isinstance(smiles, str):
        raise HTTPException(status_code=400, detail=f"Input must be a single SMILES string. Got type: `{type(smiles)}` with SMILES = `{smiles}`")

    try:
        return RDKitService().canonicalize_smiles(smiles)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/pdb/validate", tags=['Shared'])
async def validate_pdb(
    pdb_file: UploadFile = File(...),
):
    pdb_file_content = await pdb_file.read()
    try:
        return is_valid_pdb_file(pdb_file_content)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))