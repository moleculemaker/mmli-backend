import requests
import xml.etree.ElementTree as ET
import time
import json

class PubChemService:
    def __init__(self) -> None:
        pass

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

                MAX_RETRY = 20
                iteration_count = 0

                while iteration_count < MAX_RETRY:
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
                    iteration_count += 1
            else:
                print(f"SOAP Request failed with status code {smile_to_cid_response.status_code}")
