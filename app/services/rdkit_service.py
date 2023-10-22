from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

class RDKitService:
    def __init__(self) -> None:
        pass

    def renderSVGFromSMILE(self, smileString):
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

    def getFingerprint(self, smileString):
        mol = Chem.MolFromSmiles(smileString)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        return fp.ToBitString()

    def getTanimotoSimilarity(self, bitString1, bitString2):
        bitVect1 = DataStructs.CreateFromBitString(bitString1)
        bitVect2 = DataStructs.CreateFromBitString(bitString2)
        similarity = DataStructs.TanimotoSimilarity(bitVect1, bitVect2)
        return similarity