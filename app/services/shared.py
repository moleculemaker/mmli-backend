from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFMCS
from typing import Callable, List, Literal, Optional, Tuple, TypedDict

class FragmentMatchResult(TypedDict):
    mol: Chem.Mol
    matches: list[list[int]]
    # svg: str
    error: Optional[str]
    
class MCSResult(TypedDict):
    mol: Chem.Mol
    score: float
    error: Optional[str]
    
class TanimotoResult(TypedDict):
    mol: Chem.Mol
    score: float
    error: Optional[str]
    
ChemInfoAlgorithms = Literal["tanimoto", "fragment", "mcs"]
class ChemInfoResult(TypedDict):
    tanimoto: Optional[list[TanimotoResult]]
    fragment: Optional[list[FragmentMatchResult]]
    mcs: Optional[list[MCSResult]]
    
SmilesMolDict = dict[str, Chem.Mol]
MolSmilesDict = dict[Chem.Mol, str]

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
    

def convert_to_rdkit_mols(smiles_list: list[str]) -> Tuple[SmilesMolDict, MolSmilesDict]:
    """
    Convert a list of SMILES strings to a dictionary of SMILES strings and their corresponding RDKit molecules.
    
    Args:
        smiles_list (list[str]): A list of SMILES strings.
        
    Returns:
        Tuple[SmilesMolDict, MolSmilesDict]: A tuple of a dictionary of SMILES strings and their corresponding RDKit molecules, 
            and a dictionary of RDKit molecules and their corresponding SMILES strings.
            Returns None for any invalid SMILES strings in the list.
    """
    mols = { smiles: Chem.MolFromSmiles(smiles) for smiles in smiles_list }
    
    return (
        mols, 
        { mol: smiles for smiles, mol in mols.items() if mol is not None }
    )


def rdkit_get_cheminfo_algorithms_results(
    algorithms: list[ChemInfoAlgorithms],
    query_mol: Chem.Mol,
    mols: list[Chem.Mol]
) -> ChemInfoResult:
    """
    Get the results of the cheminfo algorithms for a query molecule against a list of molecules.
    
    Args:
        algorithms (list[ChemInfoAlgorithms]): The algorithms to run.
        query_mol (Chem.Mol): The query molecule.
        mols (list[Chem.Mol]): The list of molecules to compare against.
        
    Returns:
        ChemInfoResult: A dictionary of the algorithms and their results.
    """
    results = {}
    if "tanimoto" in algorithms:
        results["tanimoto"] = rdkit_get_tanimoto_similarity_score_against(query_mol, mols)
        
    if "fragment" in algorithms:
        results["fragment"] = rdkit_get_fragment_matches_against(query_mol, mols)
        
    if "mcs" in algorithms:
        results["mcs"] = rdkit_get_mcs_score_against(query_mol, mols)
        
    return results

def rdkit_get_tanimoto_similarity_score_against(
    query_mol: Chem.Mol, 
    mols: List[Chem.Mol]
) -> List[TanimotoResult]:
    """
    Get a similarity score for a query molecule against a list of molecules.
    
    Args:
        query_mol (Chem.Mol): The query molecule.
        mols (List[Chem.Mol]): The list of molecules to sort.
        
    Returns:
        List[TanimotoResult]: A list of TanimotoResult objects.
            Returns error string for any molecules failed to process.
            
    Raises:
        ValueError: If the query SMILES string is invalid
    """
    
    query_fp = Chem.RDKFingerprint(query_mol)
    
    ret_val = []
    for mol in mols:
        try:
            fp = Chem.RDKFingerprint(mol)
            similarity = DataStructs.FingerprintSimilarity(query_fp, fp)
            ret_val.append(TanimotoResult(mol=mol, score=similarity))
        except Exception as e:
            ret_val.append(TanimotoResult(mol=mol, score=-1, error=str(e)))
        
    return ret_val

    
def rdkit_get_fragment_matches_against(
    query_mol: Chem.Mol, 
    mols: list[Chem.Mol]
) -> list[FragmentMatchResult]:
    """
    Get the fragments of the query molecule in the list of molecules.
    
    Args:
        query_mol (Chem.Mol): The query molecule.
        mols (list[Chem.Mol]): The list of molecules to search for fragments.
        
    Returns:
        list[FragmentMatchResult]: A list of FragmentMatchResult objects.
            Returns error string for any molecules failed to process.
    """
    def rdkit_draw_mol_svg(mol, highlightAtoms=None, **kwargs) -> str:
        drawer = Draw.MolDraw2DSVG(300, 150)
        drawer.DrawMolecule(mol, highlightAtoms=highlightAtoms, **kwargs)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    
    ret_val = []
    for mol in mols:
        try:
            if mol.HasSubstructMatch(query_mol):
                matches = mol.GetSubstructMatches(query_mol)
                # highlight_atoms = [atom_idx for match in matches for atom_idx in match] 
                ret_val.append(FragmentMatchResult(
                    mol=mol, 
                    matches=[list(m) for m in matches],
                    # svg=rdkit_draw_mol_svg(mol, highlightAtoms=highlight_atoms)
                ))
        except Exception as e:
            ret_val.append(FragmentMatchResult(
                mol=mol, 
                matches=[], 
                # svg="", 
                error=str(e)
            ))
    return ret_val


def rdkit_get_mcs_score_against(
    query_mol: Chem.Mol, 
    mols: list[Chem.Mol]
) -> List[MCSResult]:
    """
    Get the Maximum Common Substructure (MCS) score for a query molecule against a list of molecules.
    
    Args:
        query_mol (Chem.Mol): The query molecule.
        mols (list[Chem.Mol]): The list of molecules to compare against.
        
    Returns:
        List[MCSResult]: A list of MCSResult objects.
            Returns error string for any molecules failed to process.
            The score is calculated as: MCS_atoms / (query_atoms + target_atoms - MCS_atoms)
    """
    query_atoms = query_mol.GetNumAtoms()
    ret_val = []
    
    for mol in mols:
        try:
            # Find MCS between query and current molecule
            mcs = rdFMCS.FindMCS(
                [query_mol, mol],
                completeRingsOnly=True,
                bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                timeout=1
            )
                
            # Calculate MCS score
            target_atoms = mol.GetNumAtoms()
            mcs_score = mcs.numAtoms / (query_atoms + target_atoms - mcs.numAtoms)
            ret_val.append(MCSResult(mol=mol, score=mcs_score))
            
        except Exception as e:
            ret_val.append(MCSResult(mol=mol, score=-1, error=str(e)))
            
    return ret_val



# Testing
if __name__ == "__main__":
    query_smiles = "COC1=C(C=CC(=C1)CCN)O"
    target_smiles = ["COc1ccc2c(c1OC)C(O)O[C@@H]2[C@H]1c2c(cc3c(c2OC)OCO3)CCN1C"]
    
    query_mol = Chem.MolFromSmiles(query_smiles)
    target_mols = [Chem.MolFromSmiles(smiles) for smiles in target_smiles]
    
    print(rdkit_get_tanimoto_similarity_score_against(query_mol, target_mols))
    print(rdkit_get_fragment_matches_against(query_mol, target_mols))
    print(rdkit_get_mcs_score_against(query_mol, target_mols))