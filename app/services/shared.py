from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from typing import Callable, Optional

def smiles_to_svg(smiles_string: str, 
                  width: int = 300, 
                  height: int = 150, 
                  beforeDraw: Optional[Callable[[Draw.MolDraw2DSVG], None]] = None,
                  **kwargs) -> str:
    """
    Convert a SMILES string to an SVG image.
    
    Parameters:
        smiles_string (str): SMILES string representing a molecule or a reaction
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
        if ">" in smiles_string:
            # Process as a reaction
            reaction = AllChem.ReactionFromSmarts(smiles_string, useSmiles=True)
            
            # Compute 2D coordinates for each molecule in the reaction
            for mol in reaction.GetReactants():
                rdDepictor.Compute2DCoords(mol)
                
            for mol in reaction.GetProducts():
                rdDepictor.Compute2DCoords(mol)
            
            drawer.DrawReaction(reaction, **kwargs)
        else:
            # Process as a molecule
            mol = Chem.MolFromSmiles(smiles_string)
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles_string}")
            rdDepictor.Compute2DCoords(mol)
            drawer.DrawMolecule(mol, **kwargs)

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace('svg:', '')  # Remove 'svg:' namespace for compatibility
    
    except Exception as e:
        raise Exception(f"Error processing SMILES string: {smiles_string}") from e