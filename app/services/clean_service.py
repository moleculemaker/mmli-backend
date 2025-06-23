import base64
import csv
import json
import re
from io import StringIO
import os

from fastapi import HTTPException
from sqlmodel.ext.asyncio.session import AsyncSession

from config import get_logger
from services.minio_service import MinIOService

log = get_logger(__name__)


class CleanService:

    @staticmethod
    def getProteinSeqFromDNA(DNA_sequence):
        # define a codon table as a dictionary
        codon_table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        }

        nucleotide_seq = DNA_sequence

        # assert that the length of the nucleotide sequence is a multiple of 3
        if len(nucleotide_seq) % 3 != 0:
            log.error(f'''Length of nucleotide sequence is not a multiple of 3''')
            return ''
        # split the nucleotide sequence into codons
        codons = [nucleotide_seq[i:i + 3] for i in range(0, len(nucleotide_seq), 3)]
        # translate each codon into its corresponding amino acid
        amino_acids = [codon_table[codon] for codon in codons]
        # join the resulting amino acid sequence into a string
        amino_acid_seq = "".join(amino_acids)
        # print the resulting amino acid sequence
        return amino_acid_seq


    @staticmethod
    def build_clean_job_command(job_id, job_info):
        job_config = ""

        sequence_count = 0
        # Convert JSON to FASTA format
        for record in job_info['input_fasta']:
            # Check for valid amino acid characters
            if re.match('^[ACDEFGHIKLMNPQRSTVWY]+$', record["sequence"]):
                job_config += ">{}\n{}\n".format(record["header"], record["sequence"])
            elif re.match('^[ACGTURYSWKMBDHVN]*$', record["DNA_sequence"]):
                protein_sequence = CleanService.getProteinSeqFromDNA(record["DNA_sequence"])
                if protein_sequence == '':
                    # Invalid DNA Sequence, return 400
                    raise HTTPException(status_code=400, detail='400: Bad Request - Invalid DNA Sequence')

                else:
                    job_config += ">{}\n{}\n".format(record["header"], protein_sequence)
            else:
                # Invalid FASTA Sequence, return 400
                raise HTTPException(status_code=400, detail='400: Bad Request - Invalid FASTA Protein Sequence')

            sequence_count += 1

        if sequence_count > 20:
            # Invalid FASTA Sequence, return 400
            raise HTTPException(status_code=400, detail='400: Bad Request - CLEAN allows only for a maximum of 20 FASTA Sequences.')

        encoded_data = base64.b64encode(job_config.encode('utf-8')).decode('utf-8')
        command = f'''echo {encoded_data} | base64 -d > ./data/inputs/{job_id}.fasta && (((python CLEAN_infer_fasta.py --fasta_data {job_id} 2>&1 | tee $JOB_OUTPUT_DIR/log) && ls -al ./results/inputs/; mv ./results/inputs/{job_id}_maxsep.csv $JOB_OUTPUT_DIR/{job_id}_maxsep.csv) || (touch $JOB_OUTPUT_DIR/error && false))'''

        return command


    @staticmethod
    async def cleanDBMepEsmResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """
        Outputs stored in Minio: /{job_id}/out/*  Bucket name: oed-*
        """
        folder_path = f"/{job_id}/out/"
        objects = service.list_files(bucket_name, folder_path)

        # Iterate over folder and add all contents to a dictionary
        content = {}
        for obj in objects:
            file_name = os.path.basename(obj.object_name).split('/')[-1]
            if file_name.endswith('.json'):
                content[file_name] = json.loads(service.get_file(bucket_name=bucket_name, object_name=obj.object_name))
            elif file_name.endswith('.csv'):
                content[file_name] = service.get_file(bucket_name=bucket_name, object_name=obj.object_name)
            else:
                log.warning(f'Skipping unrecognized file extension: ' + str(file_name))

        # Return the dictionary if it has contents
        if not content:
            raise HTTPException(status_code=404, detail=f"No output files were found")


    @staticmethod
    async def cleanResultPostProcess(bucket_name, job_id, service, db):
        result_path = f"{job_id}/out/{job_id}_maxsep.csv"
        csv_file = service.get_file(bucket_name, result_path)
        if csv_file is None:
            # Can't find result file, return 404
            log.error(f'CLEAN result file not found: {result_path}')
            raise HTTPException(status_code=404, detail='404: Not Found - ' + result_path)

        #log.debug(f'Returning CLEAN result file: {fullPath}')
        #input_str = ''
        #with open(fullPath, 'r') as f:
        #    input_str = f.read().strip()

        csv_lines = StringIO(csv_file.decode())
        lines = csv.reader(csv_lines, quotechar='"', delimiter=',', lineterminator='\r\n', skipinitialspace=True)
        # create a list to store the results
        results = []

        # process each line
        for line in lines:
            # split the line into the sequence header and the EC numbers/scores
            seq_header = line[0]
            ec_scores_str = line[1]

            # create a dictionary to store the sequence header and the EC numbers/scores
            result = {
                "sequence": seq_header,
                "result": []
            }

            # split the EC numbers/scores string into individual EC number/score pairs
            ec_scores_pairs = ec_scores_str.split(',')

            # process each EC number/score pair
            for ec_score_pair in ec_scores_pairs:
                # split the EC number/score pair into the EC number and the score
                ec_number, score = ec_score_pair.split('/')

                # create a dictionary to store the EC number and the score
                ec_score_dict = {
                    "ecNumber": ec_number,
                    "score": score
                }

                # add the EC number/score dictionary to the result list
                result["result"].append(ec_score_dict)

            # add the result dictionary to the results list
            results.append(result)

        # convert the results list to JSON and print it
        # json_str = json.dumps(results, indent=4)

        return json.dumps(results)
