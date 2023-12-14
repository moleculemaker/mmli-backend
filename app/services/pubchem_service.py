import requests
from sqlmodel.ext.asyncio.session import AsyncSession
from models.sqlmodel.db import get_session
from fastapi import Depends
from sqlmodel import select
from models.sqlmodel.models import MoleculeCacheEntry
import xml.etree.ElementTree as ET
import time
import json

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
    

    async def parseCidsAndGetData(self, cidFileContent, smile_list):
        # cidFileContent = str(cidFileContent)
        # print('=== CID File Content === ', cidFileContent, '===')
        molecules = cidFileContent.split('\n')
        # Remove last entry which is just a blank
        molecules.pop()
        cids_list = []
        smile_cid_dict = {}
        for molecule in molecules:
            mol_cid = molecule.split('\t')
            cids_list.append(mol_cid[1])
            smile_cid_dict[mol_cid[0]] = mol_cid[1]

        pub_chem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/MolecularFormula,MolecularWeight,IUPACName/JSON'.format(cids=",".join(cids_list))

        # Request to get data for all molecules
        response = requests.get(pub_chem_url)

        vals = []
        if response.status_code == 200:
            parsed_data = json.loads(response.text)
            properties = parsed_data['PropertyTable']['Properties']

            json_counter = 0
            for smile in smile_list:
                if smile_cid_dict[smile] == '':
                    # CID Not Available for SMILE
                    vals.append(smile)
                    vals.append('Unavailable')
                    vals.append('Unavailable')
                    vals.append('Unavailable')
                    vals.append('Unavailable')
                else:
                    vals.append(smile)
                    json_element = properties[json_counter]
                    vals.append(str(json_element["CID"]))
                    vals.append(json_element.get("MolecularFormula", "Unavailable"))
                    vals.append(json_element.get("MolecularWeight", "Unavailable"))
                    vals.append(json_element.get("IUPACName", "Unavailable"))
                    json_counter+=1
            return vals
        else:
            pass

    async def getDataForAllMolecules(self, smile_list):
            # Parse the XML data
            smile_to_cid_request_xml = '''
            <PCT-Data>
                <PCT-Data_input>
                    <PCT-InputData>
                        <PCT-InputData_query>
                            <PCT-Query>
                                <PCT-Query_type>
                                    <PCT-QueryType>
                                        <PCT-QueryType_id-exchange>
                                            <PCT-QueryIDExchange>
                                                <PCT-QueryIDExchange_input>
                                                    <PCT-QueryUids>
                                                        <PCT-QueryUids_smiles>
                                                        </PCT-QueryUids_smiles>
                                                    </PCT-QueryUids>
                                                </PCT-QueryIDExchange_input>
                                                <PCT-QueryIDExchange_operation-type value="same"/>
                                                <PCT-QueryIDExchange_output-type value="cid"/>
                                                <PCT-QueryIDExchange_output-method value="file-pair"/>
                                                <PCT-QueryIDExchange_compression value="none"/>
                                            </PCT-QueryIDExchange>
                                        </PCT-QueryType_id-exchange>
                                    </PCT-QueryType>
                                </PCT-Query_type>
                            </PCT-Query>
                        </PCT-InputData_query>
                    </PCT-InputData>
                </PCT-Data_input>
            </PCT-Data>
            '''

            job_status_xml = '''
            <PCT-Data>
            <PCT-Data_input>
                <PCT-InputData>
                <PCT-InputData_request>
                    <PCT-Request>
                    <PCT-Request_type value="status"/>
                    </PCT-Request>
                </PCT-InputData_request>
                </PCT-InputData>
            </PCT-Data_input>
            </PCT-Data>
            '''

            smile_to_cid_request_root = ET.fromstring(smile_to_cid_request_xml)
            for smile in smile_list:
                appending_smile = ET.Element("PCT-QueryUids_smiles_E")
                appending_smile.text = smile
                smile_to_cid_request_root.find(".//PCT-QueryUids_smiles").append(appending_smile)

            pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi'
            headers = {'Content-Type': 'application/xml'}
            # Request to start SMILE to CID job
            smile_to_cid_response = requests.post(pubchem_url, data=ET.tostring(smile_to_cid_request_root).decode(), headers=headers)
            if smile_to_cid_response.status_code == 200:
                smile_to_cid_response_root = ET.fromstring(smile_to_cid_response.text)
                waiting_reqid = smile_to_cid_response_root.find(".//PCT-Waiting_reqid").text

                # Creating new XML for job status fetching
                job_status_root = ET.fromstring(job_status_xml)
                request_element = job_status_root.find(".//PCT-Request")
                new_reqid = ET.Element("PCT-Request_reqid")
                new_reqid.text = waiting_reqid
                request_element.insert(0, new_reqid)

                # Convert the modified XML tree back to a string
                job_status_request_xml = ET.tostring(job_status_root).decode()

                while True:
                    # Short wait to let the job finish
                    JOB_STATUS_POLLING_PERIOD = 3
                    time.sleep(JOB_STATUS_POLLING_PERIOD)

                    # Request to check the status of the Smile to CID job
                    job_status_response = requests.post(pubchem_url, data=job_status_request_xml, headers=headers)
                    job_status_response_root = ET.fromstring(job_status_response.text)

                    status_value = job_status_response_root.find(".//PCT-Status").attrib['value']

                    if status_value == "success":
                        download_url = job_status_response_root.find(".//PCT-Download-URL_url").text
                        download_url = download_url.replace('ftp', 'https')
                        print(download_url)

                        # Request to get the file with Smile to CID conversion
                        cid_file_response = requests.get(download_url)

                        if cid_file_response.status_code == 200:
                            return await self.parseCidsAndGetData(cid_file_response.text, smile_list)
                        else:
                            pass
            else:
                print(f"SOAP Request failed with status code {smile_to_cid_response.status_code}")
