from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from typing import Callable, Optional

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