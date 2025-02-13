import argparse
import os
import re
import pyterrier as pt
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
from safetensors import safe_open
from rdkit.DataStructs.cDataStructs import CreateFromBitString

index_store_path = os.path.join(os.path.dirname(__file__), 'index_store')

########################################################
# Tokenizer for BM25
########################################################
def split_string_with_regex(text):
    KEYWORDS = ['mono','di','tri','tetra','penta','hexa','hepta','octa','nona',
            'deca','oxo','methyl','hydroxy','benzene','oxy','chloro','cyclo',
            'amino','bromo','hydro','fluoro', 'methoxy', 'ethoxy', 'phenoxy',
            'methane','cyano','amido','ethene','phospho','amide','butane',
            'carbono','hydro','sulfane','butane','sulfino', 'ethano', 'methano', 'din',
            'iodo','ethane','ethyne','bi','tri', 'imino','nitro','butan','idene','sulfo', 'propano',
            'carbon','propane','ethen','acetaldehyde','benzo','oxa','nitroso',
            'hydra','iso']
    # Create a regex pattern for the keywords
    keyword_pattern = '|'.join(sorted(KEYWORDS, key=len, reverse=True))
    
    text = text.lower()
    # Complete regex pattern
    pattern = rf'''
        {keyword_pattern}            # Match any keyword
        | \d+                        # Match standalone numbers
        | [a-zA-Z]+                  # Match alphabetic words
        | [^\s,\-\[\]\(\)\{{\}}]+    # Match any word not space, comma, hyphen, or brackets
        | (?<=\d)[a-zA-Z]+           # Match letters after a digit in the same word
    '''
    
    # Use regex to find all matches
    result = re.findall(pattern, text, re.VERBOSE)
    result = ' '.join(result)
    return result


########################################################
# SMI INDEXER
########################################################
class RDKitSubStructureSearch(pt.Transformer):
    def __init__(self, index):
        super().__init__()
        self.index = index
    def transform(self, queries, topk=1000):
        results = []
        
        meta = self.index.getMetaIndex()
        lexicon = self.index.getLexicon()

        for _,row in tqdm(queries.iterrows(), total=len(queries), desc='Processing query: '):
            qid = row['qid']
            query = row['query']
            query = Chem.MolFromSmiles(query)
            
            
            for doc_id in range(self.index.getCollectionStatistics().getNumberOfDocuments()):
                docno = meta.getItem('docno', doc_id)
                docsmi = meta.getItem('text', doc_id)
                docsmip = Chem.MolFromSmiles(docsmi)
                # Try-Catch to Ignore Invalid SMILES in Index
                try:
                    s_match = docsmip.HasSubstructMatch(query)
                    if s_match:
                        results.append({'qid':qid, 'docno':docno, 'smiles':docsmi, 'score':1.0})
                except:
                    continue
        
        result = pd.DataFrame(results)
        fil_results = pd.DataFrame()

        # Now filter the results to only topk per query
        qids = list(set(result['qid'].tolist()))
        qids = sorted(list(map(int, qids)))
        qids = list(map(str, qids))
        for qid in qids:
            qid_df = result.loc[result['qid'] == qid].iloc[:topk]
            fil_results = pd.concat([fil_results, qid_df])
        
        return fil_results, qids


def tani_search(smi_path, ten_path, query='C1COCCN1', topk=10):
    df_smi = pd.read_csv(smi_path, sep='\t')
    # Convert any ints to strings
    int_columns = df_smi.select_dtypes(include='int64').columns
    df_smi[int_columns] = df_smi[int_columns].astype(str)
    df_smi.columns = ['docno', 'text']
    
    with safe_open(ten_path, framework='pt', device='cpu') as f:
        fp_ten = f.get_tensor('fps')

    fp_ten = fp_ten.numpy()
    
    fps = []
    for fp in fp_ten:
        fps.append(CreateFromBitString(''.join(map(str, fp))))
    
    # Query Mode
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)    
    query = Chem.MolFromSmiles(query)
    q_fp = mfpgen.GetFingerprint(query)
    
    scores = DataStructs.BulkTanimotoSimilarity(q_fp, fps)
    scores = np.array(scores)
    sort_idx = np.argsort(scores)[::-1][:topk]
    scores = scores[sort_idx]
    
    # Get the SMILES for the indexes
    cand_smis = df_smi.iloc[sort_idx].reset_index(drop=True)

    return cand_smis


def rdkit_search(smi_path, query='C1COC1'):
    
    pt.java.set_log_level('DEBUG')
    # Load the pandas dataframe
    df_smi = pd.read_csv(smi_path, sep='\t')
    # df_dsmi = pd.read_csv(dsmi_path, sep='\t')

    # Convert any ints to strings
    int_columns = df_smi.select_dtypes(include='int64').columns
    df_smi[int_columns] = df_smi[int_columns].astype(str)
    df_smi.columns = ['docno', 'text']
    # int_columns = df_dsmi.select_dtypes(include='int64').columns
    # df_dsmi[int_columns] = df_dsmi[int_columns].astype(str)
    # df_dsmi.columns = ['docno', 'text', 'IUPAC']

    # Create and save the indexes
    indexer_par = pt.IterDictIndexer(index_store_path + "/smi_index", tokeniser='identity', overwrite=True,
                                     stemmer=None, stopwords=None, verbose=True,
                                     meta={'docno':20, 'text':400},
                                     properties={'max.term.length': '150'},
                                     )
    indexer_par.setProperty("max.term.length", "149")
    indexer_par.setProperty("indexer.meta.reverse.keys", "text")
    
    # indexer_diag = pt.IterDictIndexer(index_store_path + '/smid_index', overwrite=True, verbose=True,
    #                              stemmer=None, stopwords=None, tokeniser='identity',
    #                              meta={'docno':20, 'text':500, 'IUPAC':1000})
    
    indexref_par = indexer_par.index(df_smi.to_dict(orient='records'))
    index_par = pt.IndexFactory.of(indexref_par)
    
    subsearch = RDKitSubStructureSearch(index_par)
    
    queries = pd.DataFrame([
        {'qid':"1", 'query':query}
    ])

    results = subsearch.transform(queries)
    # get_paragraphs(results)
    return results
    


def get_paragraphs(results):
    par_df = 'outputs/pd_frames/Pars.tsv'
    par_df = pd.read_csv(par_df, sep='\t')
    smin_df = 'outputs/pd_frames/Par_SMI_Names.tsv'
    smin_df = pd.read_csv(smin_df, sep='\t')

    for _,row in results.iterrows():
        smi_idx = int(row['docno'])
        smi = row['smiles']
        print(f'\n For result SMILES: {smi} \n')
        names = smin_df.loc[smin_df['IDX_SMI'] == smi_idx, 'Name'].values.tolist()
        pars = par_df.loc[par_df['IDX_SMI'] == smi_idx, 'text'].values.tolist()
        pdfs = par_df.loc[par_df['IDX_SMI'] == smi_idx, 'Doc'].values.tolist()
        pages = par_df.loc[par_df['IDX_SMI'] == smi_idx, 'Page'].values.tolist()
        
        for i in range(len(names)):
            print(f'\n Found Name: {names[i]}')
        for i in range(len(pars)):
            print(f'\n Paragraph: {pars[i]}')
            print(f'\n PDF: {pdfs[i]}')
            print(f'\n Page: {pages[i]}')
        exit()            
        

########################################################
# TEXT RETRIEVAL
########################################################
def txt_eval(par_path, par2smi_path, par2img_path, smi2name_path, smi_path,
             topk=20, mm=False, cutoff=1, txt_q='chlorine'):
    # Init Pyterrier
    if not pt.started():
        pt.init()
    # Load the pandas dataframe
    df = pd.read_csv(par_path, sep='\t')

    int_columns = df.select_dtypes(include='int64').columns
    df[int_columns] = df[int_columns].astype(str)
    # Tokenize the text
    tokenize_text(df)

    
    txt_qs = pd.DataFrame([['1', split_string_with_regex(txt_q)], 
                        ], columns=['qid', 'query'])

    # Now generate the Main Index
    indexer = pt.IterDictIndexer(index_store_path + '/pd_index', overwrite=True, verbose=True, 
                                 stemmer=None, stopwords=None, tokeniser='UTFTokeniser',
                                 meta={'docno':20, 'text':8192, 'Page':3, 'Docno':50,
                                       'x1':5, 'y1':5, 'x2':5, 'y2':5})
    # Convert index df to Dict
    indexref = indexer.index(map_document(df))
    index = pt.IndexFactory.of(indexref)

    # Perform retrieval
    br = pt.BatchRetrieve(index, wmodel='BM25', num_results=topk,
                          metadata=['docno', 'Page', 'Docno', 'text', 'x1', 'y1', 'x2', 'y2'])
    results = br.transform(txt_qs)
    
    # Load the par2smi and par2img dfs
    par2smi_df = pd.read_csv(par2smi_path, sep='\t')
    int_columns = par2smi_df.select_dtypes(include='int64').columns
    par2smi_df[int_columns] = par2smi_df[int_columns].astype(str)
    
    par2img_df = pd.read_csv(par2img_path, sep='\t')
    int_columns = par2img_df.select_dtypes(include='int64').columns
    par2img_df[int_columns] = par2img_df[int_columns].astype(str)

    smi2name_df = pd.read_csv(smi2name_path, sep='\t')
    int_columns = smi2name_df.select_dtypes(include='int64').columns
    smi2name_df[int_columns] = smi2name_df[int_columns].astype(str)

    smi_df = pd.read_csv(smi_path, sep='\t')
    int_columns = smi_df.select_dtypes(include='int64').columns
    smi_df[int_columns] = smi_df[int_columns].astype(str)

    # Transform results into the form required and get the images required
    if not mm:
        results = results[['docno', 'Docno', 'Page', 'text', 'x1', 'y1', 'x2', 'y2']]
    else:
        results = results[['docno', 'Docno', 'Page', 'text', 'x1', 'y1', 'x2', 'y2', 'score']]
    img_results = pd.DataFrame(columns=['docno', 'Page', 'Name', 'SMI', 'x1', 'y1', 'x2', 'y2'])

    for _,row in results.iterrows():
        docno = row['docno']
        # Find if there are any Image SMILES attached to the passage
        docsmis_df = par2img_df.loc[par2img_df['IDX_Par'] == docno]
        if len(docsmis_df):
            for idx, (_,row) in enumerate(docsmis_df.iterrows()):
                idx_smi = row['IDX_SMI']
                smi_name = smi2name_df.loc[smi2name_df['IDX_SMI'] == idx_smi, 'name'].tolist()[0]
                smi_val = smi_df.loc[smi_df['docno'] == idx_smi, 'SMILES'].tolist()[0]
                page = docsmis_df.iloc[idx]['page']
                x1, y1 = docsmis_df.iloc[idx]['x1'], docsmis_df.iloc[idx]['y1']
                x2, y2 = docsmis_df.iloc[idx]['x2'], docsmis_df.iloc[idx]['y2']
                img_row = {'docno':docno, 'Page': page, 'Name':smi_name, 'SMI':smi_val, 'x1':x1, 'y1':y1, 'x2':x2, 'y2':y2}
                img_results = pd.concat([img_results, pd.DataFrame([img_row])])
    return results, img_results        
    


def do_sub_search(index_smi, smi_qs, par2smis, df, topk=1000):
    # Instantiate the RDKit SubStructure Search
    sub_search = RDKitSubStructureSearch(index_smi)
    smi_results, qids = sub_search.transform(smi_qs)
    # Get the relevant paragraphs per query
    smipar_df = pd.DataFrame(columns=['qid', 'docno', 'Page', 'Docno', 'text', 'x1', 'y1', 'x2', 'y2'])
    for qid in tqdm(qids, total=len(qids), desc='Populating Pars for SMI Queries: '):
        qid_smis = smi_results.loc[smi_results['qid'] == qid]
        for _, row in qid_smis.iterrows():
            smi_id = row['docno']
            par_id = par2smis.loc[par2smis['IDX_SMI'] == smi_id, 'IDX_Par'].tolist()[0]
            par_row = df.loc[df['docno'] == par_id]
            smipar_df = pd.concat([smipar_df, pd.DataFrame([{'qid':qid, 'docid': par_row.iloc[0]['docno'], 
                                                                'docno': par_row.iloc[0]['docno'], 'Page': par_row.iloc[0]['Page'],
                                                                'Docno': par_row.iloc[0]['Docno'], 'text': par_row.iloc[0]['text'],
                                                                'x1': par_row.iloc[0]['x1'], 'y1': par_row.iloc[0]['y1'],
                                                                'x2': par_row.iloc[0]['x2'], 'y2': par_row.iloc[0]['y2'],
                                                                'score': 1.0}])])
    
    return smipar_df


def do_tani_search(df_smi, smi_qs, par2smis, df, topk=20):
    # Compile the fingerprints for for the Index SMILES:
    fps = []
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True) 
    print('Computing Morgan FPs for all Index SMILES')
    for _, row in tqdm(df_smi.iterrows(), total=len(df_smi), desc='Processing... '):
        docno = row['docno']
        docsmi = row['text']
        docmol = Chem.MolFromSmiles(docsmi)
        try:
            doc_fp = mfpgen.GetFingerprint(docmol)
        except:
            doc_fp = mfpgen.GetFingerprint(Chem.MolFromSmiles('[Na]'))
        fps.append(doc_fp)
    
    cand_smis = {}
    # Compile the top 10 matches for each query
    for _,row in tqdm(smi_qs.iterrows(), total=len(smi_qs), desc='Processing query: '):
        qid = row['qid']
        query = row['query']
        query = Chem.MolFromSmiles(query)
        q_fp = mfpgen.GetFingerprint(query)

        scores = DataStructs.BulkTanimotoSimilarity(q_fp, fps)
        scores = np.array(scores)
        sort_idx = np.argsort(scores)[::-1][:topk]
        scores = scores[sort_idx]

        # Get the SMILES for the indexes
        cand_smi_ids = sort_idx
        cand_smis[qid] = list(map(str, cand_smi_ids))

    # Now compile the paragraph list
    smipar_df = pd.DataFrame(columns=['qid', 'docno', 'Page', 'Docno', 'text', 'x1', 'y1', 'x2', 'y2', 'score'])
    for qid in tqdm(cand_smis.keys(), total=len(cand_smis.keys()), desc='Populating Pars for SMI Queries: '):
        smi_ids = cand_smis[qid]
        for i, smi_id in enumerate(smi_ids):
            par_id = par2smis.loc[par2smis['IDX_SMI'] == smi_id, 'IDX_Par'].tolist()[0]
            par_row = df.loc[df['docno'] == par_id]
            smipar_df = pd.concat([smipar_df, pd.DataFrame([{'qid':qid, 'docid': par_row.iloc[0]['docno'], 
                                                                 'docno': par_row.iloc[0]['docno'], 'Page': par_row.iloc[0]['Page'],
                                                                 'Docno': par_row.iloc[0]['Docno'], 'text': par_row.iloc[0]['text'],
                                                                 'x1': par_row.iloc[0]['x1'], 'y1': par_row.iloc[0]['y1'],
                                                                 'x2': par_row.iloc[0]['x2'], 'y2': par_row.iloc[0]['y2'],
                                                                 'score': float(topk-i)}])])
    return smipar_df


# Compile the paragraphs results that match with the specific type of SMILES search
# sub -> Sub-structure searchl; tani -> Tanimoto Sim. on Morgan FP Search
def smi_eval(par_path, smi_path, par2smi_path, par2img_path, smi2name_path,
             choice='sub', topk=20, mm=False, cutoff=1, smi_q='NCCCCC1'):
    # Init Pyterrier
    if not pt.started():
        pt.init()
    # Load the paragraphs pandas dataframe
    df = pd.read_csv(par_path, sep='\t')
    int_columns = df.select_dtypes(include='int64').columns
    df[int_columns] = df[int_columns].astype(str)

    # Load the paragraphs to SMILES dataframe
    df_p2s = pd.read_csv(par2smi_path, sep='\t')
    int_columns = df_p2s.select_dtypes(include='int64').columns
    df_p2s[int_columns] = df_p2s[int_columns].astype(str)

    # Load the SMILES dataframe
    smi_df = pd.read_csv(smi_path, sep='\t')
    int_columns = smi_df.select_dtypes(include='int64').columns
    smi_df[int_columns] = smi_df[int_columns].astype(str)
    smi_df.columns = ['docno', 'text']

    # Convert the SMILES into a PyTerrier Index
    # Create and save the indexes
    indexer_par = pt.IterDictIndexer(index_store_path + "/smi_index", tokeniser='identity', overwrite=True,
                                     stemmer=None, stopwords=None, verbose=True,
                                     meta={'docno':20, 'text':400},
                                     properties={'max.term.length': '300'},
                                     )
    indexer_par.setProperty("max.term.length", "299")
    indexer_par.setProperty("indexer.meta.reverse.keys", "text")

    indexref_par = indexer_par.index(smi_df.to_dict(orient='records'))
    index_smi = pt.IndexFactory.of(indexref_par)

    # Set the queries (Split queries if Reaction SMARTS Given or multiple Strings)
    raw_qs = re.split(r'\.{1,2}|>{1,2}|\s+', smi_q.strip())
    smi_qs = []
    for id_q, q in enumerate(raw_qs):
        smi_qs.append([str(id_q+1), q])
    
    smi_qs = pd.DataFrame(smi_qs, columns=['qid', 'query'])
    
    # Do the specific RDKit Search
    
    if choice == 'sub': # Sub-Structure Search
        topk=10
        results = do_sub_search(index_smi, smi_qs, df_p2s, df, topk=topk)
    elif choice == 'tani': # Fingerprint search
        topk=10
        results = do_tani_search(smi_df, smi_qs, df_p2s, df, topk=topk)
    
    # Change the qids to query 1 and sort by score
    results['qid'] = ['1'] * len(results)
    results_counts = results['text'].value_counts().reset_index()
    results_counts.columns = ['text', 'Repeat']

    results_fo = results.drop_duplicates(subset=['text'], keep='first')
    results = results_fo.merge(results_counts, on='text')
    results = results.sort_values(by='Repeat', ascending=False)

    
    # Delete any repetitions and keep track of the number of them
    results = results[['docno', 'Docno', 'Page', 'text', 'x1', 'y1', 'x2', 'y2', 'Repeat']]

    # Load the par2smi and par2img dfs
    par2smi_df = pd.read_csv(par2smi_path, sep='\t')
    int_columns = par2smi_df.select_dtypes(include='int64').columns
    par2smi_df[int_columns] = par2smi_df[int_columns].astype(str)
    
    par2img_df = pd.read_csv(par2img_path, sep='\t')
    int_columns = par2img_df.select_dtypes(include='int64').columns
    par2img_df[int_columns] = par2img_df[int_columns].astype(str)

    smi2name_df = pd.read_csv(smi2name_path, sep='\t')
    int_columns = smi2name_df.select_dtypes(include='int64').columns
    smi2name_df[int_columns] = smi2name_df[int_columns].astype(str)

    smi_df = pd.read_csv(smi_path, sep='\t')
    int_columns = smi_df.select_dtypes(include='int64').columns
    smi_df[int_columns] = smi_df[int_columns].astype(str)

    # Transform results into the form required and get the images required
    img_results = pd.DataFrame(columns=['docno', 'Page', 'Name', 'SMI', 'x1', 'y1', 'x2', 'y2'])

    for _,row in results.iterrows():
        docno = row['docno']
        # Find if there are any Image SMILES attached to the passage
        docsmis_df = par2img_df.loc[par2img_df['IDX_Par'] == docno]
        if len(docsmis_df):
            for idx, (_,row) in enumerate(docsmis_df.iterrows()):
                idx_smi = row['IDX_SMI']
                smi_name = smi2name_df.loc[smi2name_df['IDX_SMI'] == idx_smi, 'name'].tolist()[0]
                smi_val = smi_df.loc[smi_df['docno'] == idx_smi, 'SMILES'].tolist()[0]
                page = docsmis_df.iloc[idx]['page']
                x1, y1 = docsmis_df.iloc[idx]['x1'], docsmis_df.iloc[idx]['y1']
                x2, y2 = docsmis_df.iloc[idx]['x2'], docsmis_df.iloc[idx]['y2']
                img_row = {'docno':docno, 'Page': page, 'Name':smi_name, 'SMI':smi_val, 'x1':x1, 'y1':y1, 'x2':x2, 'y2':y2}
                img_results = pd.concat([img_results, pd.DataFrame([img_row])])
    return results, img_results


def mm_eval(par_path, smi_path, par2smi_path, par2img_path, smi2name_path, 
            txt_q='chlorine', smi_q='NCCCCC1', choice='sub', cutoff=1, topk=10, fm='simple'):
    # Perform Text Retrieval
    txt_result, _ = txt_eval(par_path, par2smi_path, par2img_path, smi2name_path, smi_path,
                          txt_q=txt_q, topk=topk, mm=True)
    # Perform SMILES Retrieval
    smi_result, _ = smi_eval(par_path, smi_path, par2smi_path, par2img_path, smi2name_path,
                             smi_q=smi_q, choice=choice)
    # Fusion Method
    qids = ['1']
    if fm == 'simple':
        txt_reranked = simple_overlap(qids, txt_result, smi_result)
    else: # Number of Overlap Priority
        txt_reranked = priority_overlap(qids, txt_result, smi_result)


    # Load the par2smi and par2img dfs
    par2smi_df = pd.read_csv(par2smi_path, sep='\t')
    int_columns = par2smi_df.select_dtypes(include='int64').columns
    par2smi_df[int_columns] = par2smi_df[int_columns].astype(str)
    
    par2img_df = pd.read_csv(par2img_path, sep='\t')
    int_columns = par2img_df.select_dtypes(include='int64').columns
    par2img_df[int_columns] = par2img_df[int_columns].astype(str)

    smi2name_df = pd.read_csv(smi2name_path, sep='\t')
    int_columns = smi2name_df.select_dtypes(include='int64').columns
    smi2name_df[int_columns] = smi2name_df[int_columns].astype(str)

    smi_df = pd.read_csv(smi_path, sep='\t')
    int_columns = smi_df.select_dtypes(include='int64').columns
    smi_df[int_columns] = smi_df[int_columns].astype(str)

    # Transform results into the form required and get the images required
    results = txt_reranked[['docno', 'Docno', 'Page', 'text', 'x1', 'y1', 'x2', 'y2']]
    img_results = pd.DataFrame(columns=['docno', 'Page', 'Name', 'SMI', 'x1', 'y1', 'x2', 'y2'])

    for _,row in results.iterrows():
        docno = row['docno']
        # Find if there are any Image SMILES attached to the passage
        docsmis_df = par2img_df.loc[par2img_df['IDX_Par'] == docno]
        if len(docsmis_df):
            for idx, (_,row) in enumerate(docsmis_df.iterrows()):
                idx_smi = row['IDX_SMI']
                smi_name = smi2name_df.loc[smi2name_df['IDX_SMI'] == idx_smi, 'name'].tolist()[0]
                smi_val = smi_df.loc[smi_df['docno'] == idx_smi, 'SMILES'].tolist()[0]
                page = docsmis_df.iloc[idx]['page']
                x1, y1 = docsmis_df.iloc[idx]['x1'], docsmis_df.iloc[idx]['y1']
                x2, y2 = docsmis_df.iloc[idx]['x2'], docsmis_df.iloc[idx]['y2']
                img_row = {'docno':docno, 'Page': page, 'Name':smi_name, 'SMI':smi_val, 'x1':x1, 'y1':y1, 'x2':x2, 'y2':y2}
                img_results = pd.concat([img_results, pd.DataFrame([img_row])])
    
    return results, img_results
    

def simple_overlap(qids, txt_result, smi_result):
    txt_reranked = pd.DataFrame()
    # Rerank txt scores based on if a match is found for the same par in SMIs
    qid_txt = txt_result
    qid_smi = smi_result
    overlap = qid_txt.loc[qid_txt['docno'].isin(qid_smi['docno'])]
    non_overlap = qid_txt.loc[~qid_txt['docno'].isin(qid_smi['docno'])]
    new_merged = pd.concat([overlap, non_overlap])
    temp_scores = np.sort(np.array(new_merged['score'].values.tolist()))[::-1]
    new_merged['score'] = temp_scores
    txt_reranked = pd.concat([txt_reranked, new_merged])

    return txt_reranked


def priority_overlap(qids, txt_result, smi_result):
    txt_reranked = pd.DataFrame()

    # Get the number of repetitions as a set
    rep_vals = sorted(list(set(smi_result['Repeat'].values.tolist())), reverse=True)
    qid_txt = txt_result
    qid_smi = smi_result

    # Merge Only the overlap values first
    for val in rep_vals:
        val_rows = smi_result.loc[smi_result['Repeat'] == val]
        overlap = qid_txt.loc[qid_txt['docno'].isin(val_rows['docno'])]
        txt_reranked = pd.concat([txt_reranked, overlap])
    
    
    # Now Merge Only all the non-overlapping values
    non_overlaps = qid_txt.loc[~qid_txt['docno'].isin(txt_reranked['docno'])]
    txt_reranked = pd.concat([txt_reranked, non_overlaps])
    temp_scores = np.sort(np.array(qid_txt['score'].values.tolist()))[::-1]
    txt_reranked['score'] = temp_scores

    return txt_reranked


def map_document(df):
    return df.to_dict(orient='records')   

def split_query_txt(text):
    return text.split('|')[0]

def split_query_smi(text):
    return text.split('|')[1]

def tokenize_text(df, col='text'):
    df[col] = df[col].apply(split_string_with_regex)
    return df

# Splits the mumtimodal queries into either txt or smiles
def split_mm_query(df, mode='txt'):
    if mode == 'txt':
        df['query'] = df['query'].apply(split_query_txt)
    else:
        df['query'] = df['query'].apply(split_query_smi)
    return df


def parse_args():
    parser = argparse.ArgumentParser(description='Argument parser for app')
    parser.add_argument('--txt_q', type=str, default='synthesis of fluorooxetan pyrimidine')
    parser.add_argument('--smi_q', type=str, default='ClC1=NC=CC=N1')
    parser.add_argument('--par_path', type=str, default='indexing_data_susuki/tsv_files/Pars.tsv')
    parser.add_argument('--smi_path', type=str, default='indexing_data_susuki/tsv_files/Par_SMI.tsv')
    parser.add_argument('--par2smi_path', type=str, default='indexing_data_susuki/tsv_files/Pars2SMI.tsv')
    parser.add_argument('--par2img_path', type=str, default='indexing_data_susuki/tsv_files/Par_Imgs.tsv')
    parser.add_argument('--smi2name_path', type=str, default='indexing_data_susuki/tsv_files/Par_SMI_Names.tsv')
    parser.add_argument('--out_dir', type=str, default='outputs')

    args = parser.parse_args()
    smi_q = args.smi_q
    txt_q = args.txt_q
    
    # Text Only Search
    if not len(smi_q):
        result, img_result = txt_eval(args.par_path, args.par2smi_path, args.par2img_path, 
                          args.smi2name_path, args.smi_path, txt_q=txt_q, topk=10)
    # SMILES Only Search
    elif not len(txt_q):
        result, img_result = smi_eval(args.par_path, args.smi_path, args.par2smi_path, 
                        args.par2img_path, args.smi2name_path, 
                        choice='tani', smi_q=smi_q)
    # Multi-modal Search
    else:
        result, img_result = mm_eval(args.par_path, args.smi_path, args.par2smi_path,
                        args.par2img_path, args.smi2name_path, txt_q=txt_q, smi_q=smi_q,
                        choice='sub', topk=10, fm='priority')
    
    # Save the results
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    
    result.to_csv(os.path.join(args.out_dir, 'passage_results.tsv'), 
                  index=False, sep='\t')
    img_result.to_csv(os.path.join(args.out_dir, 'image_results.tsv'), 
                  index=False, sep='\t')

    return result, img_result
    

if __name__ == '__main__':
    result, img_result = parse_args()