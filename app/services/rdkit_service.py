from rdkit import Chem
from rdkit.Chem import Draw

class RDKitService:
    def __init__(self) -> None:
        pass

    def renderSVGFromSMILE(smileString):
        try:
            m = Chem.MolFromSmiles(smileString)
            d = Draw.rdMolDraw2D.MolDraw2DSVG(120, 120)
            if m is None:
                raise Exception("An Error ocurred creating the molecule structure")
            d.DrawMolecule(m)
            d.FinishDrawing()
            p = d.GetDrawingText()
            return p
        except Exception as error:
            print("Error: ", error)
            return "Unavailable"