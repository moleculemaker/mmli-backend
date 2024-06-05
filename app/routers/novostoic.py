from typing import Optional, Union

from fastapi import APIRouter, Depends
from pydantic import BaseModel
from sqlalchemy import or_

from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from models.sqlmodel.db import get_session
from models.sqlmodel.models import ChemicalIdentifier

from rdkit.Chem import CanonSmiles

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import re

class ChemicalAutoCompleteResponse(BaseModel):
    name: str
    smiles: str
    inchi: str
    inchi_key: str
    metanetx_id: str
    kegg_id: str
    structure: Optional[str]

router = APIRouter()

@router.get(f"/chemical/auto-complete", 
            tags=['Novostoic'], 
            response_model=list[ChemicalAutoCompleteResponse], 
            description="Returns a list of chemicals that match the search string limited to 20 results.")
async def get_chemical_auto_complete(search: str, db: AsyncSession = Depends(get_session)):
    existing_chemicals = await db.execute(
        select(ChemicalIdentifier.name, 
               ChemicalIdentifier.smiles, 
               ChemicalIdentifier.inchi, 
               ChemicalIdentifier.inchi_key, 
               ChemicalIdentifier.metanetx_id, 
               ChemicalIdentifier.kegg_id
        ).filter(or_(
            ChemicalIdentifier.name.like(f"{search}%"),
            ChemicalIdentifier.smiles.like(f"{search}%"),
            ChemicalIdentifier.inchi.like(f"{search}%"),
            ChemicalIdentifier.inchi_key.like(f"{search}%"),
            ChemicalIdentifier.metanetx_id.like(f"{search}%"),
            ChemicalIdentifier.kegg_id.like(f"{search}%")
        )).limit(20)
    )

    return [
        {
            "name": chemical[0].lower(),
            "smiles": chemical[1],
            "inchi": chemical[2],
            "inchi_key": chemical[3],
            "metanetx_id": chemical[4],
            "kegg_id": chemical[5]
        } for chemical in existing_chemicals.all()
    ]

@router.get(f"/chemical/validate", tags=['Novostoic'], response_model=Union[ChemicalAutoCompleteResponse, None])
async def validate_chemical(search: str, db: AsyncSession = Depends(get_session)):
    try:
        smiles = CanonSmiles(search)
    except:
        smiles = search

    existing_chemicals = (await db.execute(
        select(ChemicalIdentifier.name, 
               ChemicalIdentifier.smiles, 
               ChemicalIdentifier.inchi, 
               ChemicalIdentifier.inchi_key, 
               ChemicalIdentifier.metanetx_id, 
               ChemicalIdentifier.kegg_id
        ).filter(or_(
            ChemicalIdentifier.smiles == smiles,
            ChemicalIdentifier.name == search,
            ChemicalIdentifier.inchi == search,
            ChemicalIdentifier.inchi_key == search,
            ChemicalIdentifier.metanetx_id == search,
            ChemicalIdentifier.kegg_id == search)
        )
    )).all()

    if not len(existing_chemicals):
        return 
    
    chemical = existing_chemicals[0]
    pattern = re.compile("<\?xml.*\?>")

    def draw_mol(mol, molSize=(300,150), kekulize=True):
        mc = Chem.MolFromSmiles(mol)
        if kekulize:
            try:
                Chem.Kekulize(mc)
            except:
                mc = Chem.Mol(mol.ToBinary())
        if not mc.GetNumConformers():
            Chem.rdDepictor.Compute2DCoords(mc)

        drawer = rdMolDraw2D.MolDraw2DSVG(*molSize)
        drawer.DrawMolecule(mc)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:', '')
        svg = re.sub(pattern, '', svg)
        return svg

    return {
        "name": chemical[0],
        "smiles": chemical[1],
        "inchi": chemical[2],
        "inchi_key": chemical[3],
        "metanetx_id": chemical[4],
        "kegg_id": chemical[5],
        "structure": draw_mol(chemical[1]) if chemical[1] else None
    }
