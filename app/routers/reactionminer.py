from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel

from services.reactionminer_service import ReactionMinerService

class ReactionMinerSearchResponse(BaseModel):
    status: str
    passage_results: list[dict]
    image_results: list[dict]

router = APIRouter()

@router.get(f"/reactionminer/search", 
            tags=['ReactionMiner'], 
            response_model=ReactionMinerSearchResponse,
            description="Returns a list of passages and images that match the search string.")
def search_reactionminer(smi_search: str = Query(default="", description="SMILES search string"), 
                        txt_search: str = Query(default="", description="Text search string")
    ):
    result = ReactionMinerService.search_reactionminer(smi_search, txt_search)    
    if result['status'] == 'error':
        raise HTTPException(status_code=500, detail=result['error'])
    
    return ReactionMinerSearchResponse(**result)