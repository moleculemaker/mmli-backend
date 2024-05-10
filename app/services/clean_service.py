import base64
import re

from fastapi import HTTPException

from config import get_logger

log = get_logger(__name__)


def getProteinSeqFromDNA(self, DNA_sequence):
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


def build_clean_job_command(job_id, job_info):
    job_config = ""

    sequence_count = 0
    # Convert JSON to FASTA format
    for record in job_info['input_fasta']:
        # Check for valid amino acid characters
        if re.match('^[ACDEFGHIKLMNPQRSTVWY]+$', record["sequence"]):
            job_config += ">{}\n{}\n".format(record["header"], record["sequence"])
        elif re.match('^[ACGTURYSWKMBDHVN]*$', record["DNA_sequence"]):
            protein_sequence = getProteinSeqFromDNA(record["DNA_sequence"])
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
