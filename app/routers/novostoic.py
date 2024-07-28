from typing import Optional, Union

from fastapi import APIRouter, Depends, Query
from pydantic import BaseModel
from sqlalchemy import or_, text
from sqlalchemy.orm import aliased

from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from models.sqlmodel.db import get_session
from models.sqlmodel.models import ChemicalIdentifier, ReactionEnzID, EnzIDEnzIID, EnzIIDEnzSeq

from rdkit.Chem import CanonSmiles
from services.shared import draw_chemical_svg

class ChemicalAutoCompleteResponse(BaseModel):
    name: str
    smiles: str
    inchi: str
    inchi_key: str
    metanetx_id: str
    kegg_id: str
    structure: Optional[str]

router = APIRouter()

@router.get(f"/novostoic/auto-complete", 
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

@router.get(f"/novostoic/validate", tags=['Novostoic'], response_model=Union[ChemicalAutoCompleteResponse, None])
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
               ChemicalIdentifier.kegg_id,
               ChemicalIdentifier.is_cofactor,
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
    
    
    # if there are multiple matches, return the first one that has a KEGG ID
    chemicals_with_keggid = [chemical for chemical in existing_chemicals if chemical[5]]
    chemical = chemicals_with_keggid[0] if len(chemicals_with_keggid) else existing_chemicals[0]

    return {
        "name": chemical[0],
        "smiles": chemical[1],
        "inchi": chemical[2],
        "inchi_key": chemical[3],
        "metanetx_id": chemical[4],
        "kegg_id": chemical[5],
        "is_cofactor": chemical[6] is not None,
        "structure": draw_chemical_svg(chemical[1]) if chemical[1] else None
    }
    
@router.get("/novostoic/enzseqs", tags=['Novostoic'])
async def get_metanex_reaction_enzymes(
    metanetx: str = Query(..., description="Comma-separated list of reaction IDs"),
    db: AsyncSession = Depends(get_session)
):
    # Split the comma-separated string into a list
    reaction_ids = [id.strip() for id in metanetx.split(',')]

    re = aliased(ReactionEnzID)
    ee = aliased(EnzIDEnzIID)
    es = aliased(EnzIIDEnzSeq)

    query = (
        select(re.reaction_id, re.enz_id, es.seq)
        .join(ee, re.enz_id == ee.enz_id)
        .join(es, ee.enz_iid == es.id)
        .where(re.reaction_id.in_(reaction_ids))
    )

    # Execute the query
    result = await db.execute(query, {"reaction_ids": reaction_ids})
    
    # Fetch all results
    all_results = result.all()
    
    # Group results by reaction_id
    grouped_results = {}
    for row in all_results:
        reaction_id = row.reaction_id
        if reaction_id not in grouped_results:
            grouped_results[reaction_id] = []
        grouped_results[reaction_id].append({"enz_id": row.enz_id, "seq": row.seq})
    
    return grouped_results