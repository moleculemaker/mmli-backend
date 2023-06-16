import requests

class PubChemService:
    def __init__(self) -> None:
        pass

    def queryMoleculeProperties(smile):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smile}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
        # url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/IC1=CC=CC=C1/property/MolecularFormula,MolecularWeight,IUPACName/JSON"

        response = requests.get(url)

        cid = 'Unknown'
        molecularFormula = 'Unknown'
        molecularWeight = 'Unknown'
        name = 'Unknown'
        if response.status_code == 200:
            data = response.json()
            # Extract properties
            properties = data['PropertyTable']['Properties'][0]
            cid = properties['CID']
            if cid and cid != '0':
                molecularFormula = properties['MolecularFormula']
                molecularWeight = properties['MolecularWeight']
                name = properties['IUPACName']
        else:
            print(f"Request failed with status code {response.status_code}")
        return cid, name, molecularFormula, molecularWeight

    def getAdditionalProperties(cid):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/"
        # url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/11587/JSON/"
        response = requests.get(url)

        chemicalSafety = []
        Description = 'Unknown'
        if response.status_code == 200:
            data = response.json()
            # chemicalSafety = 
            try:
                sections = data.get('Record').get('Section')

                match = next((obj for obj in sections if obj.get("TOCHeading") == "Names and Identifiers"), None)
                if match is not None:
                    nestedSections = match.get('Section')
                    nestedMatch = next((obj for obj in nestedSections if obj.get("TOCHeading") == "Record Description"), None)
                    if nestedMatch is not None:
                        information = nestedMatch.get('Information')
                        descriptionMatch = next((obj for obj in information if obj.get("Name") == "Record Description" or obj.get("Description") == "Physical Description"), None)
                        Description = descriptionMatch.get('Value').get('StringWithMarkup')[0].get('String')
            except Exception as error:
                print("Error getting molecule description: ", error)

            try:
                chemicalSafetyMatch = next((obj for obj in sections if obj.get("TOCHeading") == "Chemical Safety"), None)
                if chemicalSafetyMatch is not None:
                    markupArray = chemicalSafetyMatch.get('Information')[0].get('Value').get('StringWithMarkup')[0].get('Markup')
                    chemicalSafety = ([obj["Extra"] for obj in markupArray if "Extra" in obj])
            except Exception as error:
                print("Error getting chemical safety information: ", error)
        else:
            print(f"Request failed with status code {response.status_code} for {cid}")
        return chemicalSafety, Description