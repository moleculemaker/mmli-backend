from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem import DataStructs
from typing import Callable, Dict, List, Optional, Tuple, TypeVar, Literal

T = TypeVar('T')
ValueOrError = Tuple[Literal[0], T] | Tuple[int, str]

def draw_chemical_svg(id: str,
                  width: int = 300, 
                  height: int = 150, 
                  beforeDraw: Optional[Callable[[Draw.MolDraw2DSVG], None]] = None,
                  **kwargs) -> str:
    """
    Convert a chemical string to an SVG image.
    
    Parameters:
        id (str): string representing a molecule or a reaction, can be inchi, smiles, or reaction smarts
        width (int): Width of the SVG image
        height (int): Height of the SVG image
        beforeDraw (Callable[[Draw.MolDraw2DSVG], None], optional): Function to customize drawer before drawing
        **kwargs: Additional keyword arguments to pass to the drawing function
        
    Returns: 
        str: SVG image as a string
        
    Raises:
        Exception: If the SMILES string is invalid
    """
    try:
        drawer = Draw.MolDraw2DSVG(width, height)
        
        # Apply custom drawing options if provided
        if beforeDraw:
            beforeDraw(drawer)
        
        # Check if the SMILES string represents a reaction
        if ">" in id:
            # Process as a reaction
            reaction = AllChem.ReactionFromSmarts(id, useSmiles=True)
            
            # Compute 2D coordinates for each molecule in the reaction
            for mol in reaction.GetReactants():
                rdDepictor.Compute2DCoords(mol)
                
            for mol in reaction.GetProducts():
                rdDepictor.Compute2DCoords(mol)
            
            drawer.DrawReaction(reaction, **kwargs)
        else:
            # Process as a molecule
            if id.lower().startswith("inchi="):
                mol = Chem.MolFromInchi(id)
            else:
                mol = Chem.MolFromSmiles(id)
                
            if mol is None:
                raise ValueError(f"Invalid input: {id}")
            
            rdDepictor.Compute2DCoords(mol)
            drawer.DrawMolecule(mol, **kwargs)

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace('svg:', '')  # Remove 'svg:' namespace for compatibility
    
    except Exception as e:
        raise Exception(f"Error processing input: {id}") from e
    

def get_similarity_score(query_smiles: str, smiles_list: List[str]) -> ValueOrError[Dict[str, float | str]]:
    """
    Get a similarity score for a query SMILES string against a list of SMILES strings.
    
    Args:
        query_smiles (str): The query SMILES string.
        smiles_list (List[str]): The list of SMILES strings to sort.
        
    Returns:
        ValueOrError[Dict[str, float]]: A dictionary of SMILES strings and their similarity scores, or an exception message.
        
    Raises:
        Exception: If the query SMILES string is invalid
    """
    try:
        query_fp = Chem.RDKFingerprint(Chem.MolFromSmiles(query_smiles))
    except Exception as e:
        return (1, f"invalid query smiles: {query_smiles}")
    
    ret_val = {}
    for smiles in smiles_list:
        try:
            fp = Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))
        except Exception as e:
            ret_val[smiles] = "invalid smiles"
            continue
        
        similarity = DataStructs.FingerprintSimilarity(query_fp, fp)
        ret_val[smiles] = similarity
        
    return (0, ret_val)
    