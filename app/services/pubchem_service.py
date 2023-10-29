import requests
from sqlmodel.ext.asyncio.session import AsyncSession
from models.sqlmodel.db import get_session
from fastapi import Depends
from sqlmodel import select
from models.sqlmodel.models import MoleculeCacheEntry

class PubChemService:
    def __init__(self, db) -> None:
        self.separator = "|"
        self.db = db

    async def check_molecule_cache(self, smile):
        result = await self.db.execute(select(MoleculeCacheEntry).where(MoleculeCacheEntry.smile == smile))
        # result = await self.db.get(MoleculeCacheEntry, smile)
        molecule_entry = result.scalar()
        if molecule_entry:
            # Deserialize Chemical Safety
            chemical_safety = molecule_entry.chemical_safety.split(self.separator)
            return molecule_entry.pub_chem_id, molecule_entry.name, molecule_entry.molecular_formula, molecule_entry.molecular_weight, chemical_safety, molecule_entry.description
        else:
            raise Exception()

    async def save_molecule_to_cache(self, smile, cid, name, molecularFormula, molecularWeight, chemicalSafety, Description):
            # Serialize chemical safety array
            chemicalSafety = self.separator.join(chemicalSafety) 

            molecule_entry = MoleculeCacheEntry(
                smile = smile,
                pub_chem_id = cid,
                name = name,
                molecular_formula = molecularFormula,
                molecular_weight = molecularWeight,
                chemical_safety = chemicalSafety,
                description = Description
            )

            self.db.add(molecule_entry)
            await self.db.commit()

    async def queryMoleculeProperties(self, smile):
        try:
            return await self.check_molecule_cache(smile=smile)
        except:
            # Get from Pub Chem API
            pass

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smile}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
        # url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/IC1=CC=CC=C1/property/MolecularFormula,MolecularWeight,IUPACName/JSON"

        response = requests.get(url)

        cid = 'Unavailable'
        molecularFormula = 'Unavailable'
        molecularWeight = 'Unavailable'
        name = 'Unavailable'
        chemicalSafety = []
        Description = 'Unavailable'
        if response.status_code == 200:
            data = response.json()
            # Extract properties
            properties = data['PropertyTable']['Properties'][0]
            cid = properties['CID']
            if cid and cid != '0':
                molecularFormula = properties['MolecularFormula']
                molecularWeight = properties['MolecularWeight']
                name = properties['IUPACName']
                chemicalSafety, Description = await self.getAdditionalProperties(cid)
        else:
            print(f"Request failed with status code {response.status_code}")

        # Save to Cache
        try: 
            await self.save_molecule_to_cache(smile, cid, name, molecularFormula, molecularWeight, chemicalSafety, Description)
        except:
            print("Unable to save molecule to molecule cache")

        return cid, name, molecularFormula, molecularWeight, chemicalSafety, Description

    async def getAdditionalProperties(self, cid):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/"
        # url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/11587/JSON/"
        response = requests.get(url)

        chemicalSafety = []
        Description = 'Unavailable'
        if response.status_code == 200:
            data = response.json()
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