from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from typing import Callable, Optional
import requests
from fastapi import HTTPException

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
        
def get_iupac_name(smiles: str) -> str:
    try:
        escaped_smiles = smiles.replace('#', '%23') # "Triple bonds in SMILES strings represented by ‘#’ have to be URL-escaped as ‘%23’"
        url = f"https://cactus.nci.nih.gov/chemical/structure/{escaped_smiles}/iupac_name"
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        else:
            return None
    except Exception as e:
        return None
    
    
def is_valid_pdb_file(file_content) -> bool:
    try:
        # RDKit's MolFromPDBBlock expects the PDB block as text; callers may pass
        # raw bytes (e.g. from an upload or MinIO), so decode to str first.
        if isinstance(file_content, (bytes, bytearray)):
            file_content = file_content.decode('utf-8', errors='ignore')
        mol = Chem.MolFromPDBBlock(file_content)
        return mol is not None
    except Exception:
        return False


# Fallback when a tool's config omits `maxResidues`. Kept as a named constant so the
# bound is greppable from both call sites; the live value comes from config.
DEFAULT_MAX_STRUCTURE_RESIDUES = 700


def _residue_count(lines: list) -> int:
    """Length of one record's body: whitespace removed, trailing stop codon dropped."""
    seq = ''.join(''.join(lines).split())
    return len(seq.rstrip('*'))


def count_fasta_residues(fasta: str) -> int:
    """Residue count of the longest record in a FASTA string.

    Longest record, not the total: the constraint being guarded is peak GPU memory
    during a fold, which is driven by the largest single chain, not by how many chains
    were submitted.

    Tolerant of the shapes that actually arrive -- a bare sequence with no header,
    wrapped lines, blank lines, a trailing '*' -- because this runs before any
    tool-specific parsing and must not reject input the tool would have accepted.
    """
    longest = 0
    current = []
    for line in (fasta or '').splitlines():
        if line.startswith('>'):
            longest = max(longest, _residue_count(current))
            current = []
        else:
            current.append(line)
    return max(longest, _residue_count(current))


def validate_fasta_length(fasta: str, max_residues: int, tool_name: str) -> None:
    """Reject a FASTA whose longest record exceeds `max_residues`.

    Enforced server-side because the browser is not the only client: the frontend warns
    and omits the structure job above this bound, but a direct API call would otherwise
    put a fold on the shared GPU that cannot fit. Exceeding it is not a server fault and
    not a transient condition, so it is a 422 rather than a 400 or a retry.
    """
    residues = count_fasta_residues(fasta)
    if residues > max_residues:
        raise HTTPException(
            status_code=422,
            detail=(f'{tool_name}: sequence is {residues} residues, which exceeds the '
                    f'{max_residues}-residue limit for structure prediction'),
        )
