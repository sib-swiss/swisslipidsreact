import pandas as pd
import random
import os

from pyrheadb.RheaDB import RheaDB
from .SwissLipids import SwissLipids
from .RheaToSwisslipidsDf import RheaToSwisslipidsDf
import networkx as nx
import importlib.resources

from platformdirs import user_cache_dir
import re
from .FA_lists import positions, get_FA_list, FAS_15, FAS_85, FAS_79, FOH_15, PAL_C16, PALOH_C16, PAL_C16_OCT_C18, SPHINGO_23

class MasterIdAnalysis:
    
    def __init__(self, output_dir='', timestamp='now'):
        # get a safe cache directory
        self.cache_dir = user_cache_dir("swisslipidsreact")
        os.makedirs(self.cache_dir, exist_ok=True)
        self.output_dir = output_dir
        self.timestamp=timestamp

    def get_lipid_class_graph(self):
        # Build a directed graph of lipid class relationships
        # Split multiple lipid classes into separate rows
        df_expanded = self.swisslipids.assign(
            **{'Lipid class*': self.swisslipids['Lipid class*'].str.split('|')}
        ).explode('Lipid class*')

        # Optional: strip whitespace
        df_expanded['Lipid class*'] = df_expanded['Lipid class*'].str.strip()

        # Now build the graph
        self.G_lipid_class = nx.from_pandas_edgelist(
            df_expanded, source='Lipid class*', target='Lipid ID', create_using=nx.DiGraph()
        )

    # ---------- Lipid Class Graph Analysis ----------
    def get_lipid_to_descendant_df(self, parent_SLMs):

        # Find descendants for each lipid in the Rhea-SwissLipids merged set
        list_id_descendant = []
        for lipid_id in set(parent_SLMs):
            descendants = nx.descendants(self.G_lipid_class, lipid_id)
            list_id_descendant.extend((lipid_id, i) for i in descendants)

        # Create dataframe of isomeric subspecies relationships
        lipid_to_descendant_df = pd.DataFrame(list_id_descendant, columns=[
            'rhea_lipid_id', 'isomeric_subspecies_descendant_lipid_id'
        ])

        lipid_to_descendant_df = lipid_to_descendant_df.merge(
            self.df_isomeric_subspecies, left_on='isomeric_subspecies_descendant_lipid_id',
            right_on='Lipid ID', how='inner'
        )
        return lipid_to_descendant_df
    
    def get_only_chebi_sl_df(self):
        """
        Extracts and processes CHEBI identifiers from the SwissLipids dataframe.
        Returns a dataframe mapping 'Lipid ID' to individual numeric ChEBI IDs.
        """
        df_chebi = self.swisslipids[['Lipid ID', 'CHEBI', 'Level']].copy()
        df_chebi = df_chebi[df_chebi['Level']!='Isomeric subspecies']
        df_chebi.drop(columns=['Level'], inplace=True)
        df_chebi.dropna(subset=['CHEBI'], inplace=True)
        df_chebi['CHEBI'] = df_chebi['CHEBI'].astype(str).str.split('|')
        df_chebi = df_chebi.explode('CHEBI')
        df_chebi.dropna(subset=['CHEBI'], inplace=True)
        df_chebi['chebi_id'] = df_chebi['CHEBI'].apply(lambda x: int(float(x.replace('CHEBI:', ''))) if isinstance(x, str) else int(x))
        self.df_chebi = df_chebi
        return df_chebi
   
    def SLMs_from_CHEBIs(self, list_of_chebi_ids):

        return self.df_chebi[self.df_chebi['chebi_id'].isin(list_of_chebi_ids)]['Lipid ID'].tolist()
    
    def get_isomeric_subspecies_table(self):
        # Filter isomeric subspecies
        isomeric_subspecies = self.swisslipids[self.swisslipids['Level'] == 'Isomeric subspecies']
        self.df_isomeric_subspecies = isomeric_subspecies.copy()


    def run_master_id_analysis(self, results_overview_path='path/to/results.tsv', curated_fa_list_run=True, output_dir=None, no_curated_list_restrictions=True):
        results_overview = pd.read_csv(results_overview_path, sep='\t')
        # determine the base directory
        if output_dir is None:
            output_dir = os.getcwd()
        else:
            os.makedirs(output_dir, exist_ok=True)

        # ---------- Load and Process RheaDB Data ----------

        # Top summary lines
        summary_lines = []

        # First read SwissLipids
        sl = SwissLipids(output_dir=output_dir)
        sl.read_swisslipids_from_file()

        print('Unique SLMs', len(sl.swisslipids['Lipid ID'].unique()))
        print('Unique SLMs isomeric subspecies', len(sl.swisslipids[sl.swisslipids['Level'] == 'Isomeric subspecies']['Lipid ID'].unique()))

        r2sl=RheaToSwisslipidsDf()

        # Map SwissLipids to ChEBI, and merge with Rhea data
        df_swiss_lipids_chebi = sl.get_only_chebi_sl_df()

        summary_lines.append(f'# swiss lipids * chebi\t{len(df_swiss_lipids_chebi)}')
        print('Unique ChEBI IDs in SwissLipids:', len(df_swiss_lipids_chebi['CHEBI'].unique()))
        print('Unique SL IDs with ChEBI:', len(df_swiss_lipids_chebi['Lipid ID'].unique()))

        rdb = RheaDB()
        rheadf = rdb.rhea_reaction_long_format_smiles_chebi
        rhea_reactions = rdb.df_reactions

        rhea_non_residue_reactions = rhea_reactions[rhea_reactions['residue_rxn_flag']==False]
        rhea_class_reactions = rhea_reactions[rhea_reactions['Web-RInChIKey'].isna()]
        rhea_class_reactions = rhea_class_reactions[rhea_class_reactions['class_reaction_flag']==True]
        rheadf = rheadf[rheadf['MASTER_ID'].isin(rhea_non_residue_reactions['MASTER_ID'])]
        rheadf = rheadf[rheadf['MASTER_ID'].isin(rhea_class_reactions['MASTER_ID'])]
        rhea_reactions = rhea_reactions.copy()
        rhea_reactions['1_initial_long_df'] = rhea_reactions['MASTER_ID'].isin(rheadf['MASTER_ID'])

        #### SWISS LIPIDS CORRECTIONS ####

        # In swisslipids 58342 (ChEBI) corresponds to Fatty acyl-CoAs, 
        # however, in ChEBI it is acyl-CoA, and in swisslipids CHEBI:77636 does not exist
        # Yes, it is confusing, thing is fatty acyl is a child of acyl
        # But acyl is not always a lipid
        # However acyl was used since the beginning of times to describe the fatty acyl reactions
        # ->
        # I will change all CHEBI:77636 to CHEBI:58342 in the input to address this
        rheadf.loc[rheadf['chebiid'] == 'CHEBI:77636', 'chebiid'] = 'CHEBI:58342'

        # Replacing 2-monolysocardiolipin(2−) and 1-monolysocardiolipin(2−) with monolysocardiolipin(2−), since it is the enumerated version in SwissLipids
        rheadf.loc[rheadf['chebiid'] == 'CHEBI:64743', 'chebiid'] = 'CHEBI:167057'
        rheadf.loc[rheadf['chebiid'] == 'CHEBI:65092', 'chebiid'] = 'CHEBI:167057'

        # Replacing 2,3-diacyl-sn-glycerol with 1,2-diglyceride
        rheadf.loc[rheadf['chebiid'] == 'CHEBI:75524', 'chebiid'] = 'CHEBI:49172'

        # changing carboxilate to fatty acid for RHEA:34359 and RHEA:36195
        rheadf.loc[(rheadf['chebiid'] == 'CHEBI:29067') & (rheadf['MASTER_ID'] == 34359), 'chebiid'] = 'CHEBI:28868'
        rheadf.loc[(rheadf['chebiid'] == 'CHEBI:29067') & (rheadf['MASTER_ID'] == 36195), 'chebiid'] = 'CHEBI:28868'

        # changing aliphatic alcool to fatty alcohol for RHEA:59388
        rheadf.loc[(rheadf['chebiid'] == 'CHEBI:2571') & (rheadf['MASTER_ID'] == 59388), 'chebiid'] = 'CHEBI:24026'

        # changing a long chain fatty alcohol to fatty alcohol for RHEA:38443 since info about the chain length will come from the wax ester
        rheadf.loc[(rheadf['chebiid'] == 'CHEBI:17135') & (rheadf['MASTER_ID'] == 38443), 'chebiid'] = 'CHEBI:24026'

        # changing sterol to cholesterol since we only put cholesterol esters into the enumeration strategy
        # have to hardcode the SMILES for cholesterol (copied from SwissLipids)
        rheadf.loc[rheadf['chebiid'] == 'CHEBI:15889', 'smiles'] = 'CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C'
        rheadf.loc[rheadf['chebiid'] == 'CHEBI:15889', 'chebiid'] = 'CHEBI:16113'

        # Remove extract numeric ChEBI IDs

        def split_chebi(chebi):
            """
            Splits a ChEBI identifier (e.g., "CHEBI:1234") into prefix and numeric ID.
            Returns (prefix, numeric ID) or (None, None) on failure.
            """
            try:
                prefix, num = chebi.split(':')
                return pd.Series([prefix, int(num)])
            except Exception:
                return pd.Series([None, None])
        
        rheadf[['class', 'chebi_id']] = rheadf['chebiid'].apply(split_chebi)
        rheadf.dropna(subset=['chebi_id'], inplace=True)
        rheadf.loc[:, 'chebi_id'] = rheadf['chebi_id'].astype(int)

        # Remove polymer reactions
        df_rhea_polymer_reactions = rheadf[rheadf['class'] == 'POLYMER']
        df_rhea_chebi = rheadf[~rheadf['MASTER_ID'].isin(df_rhea_polymer_reactions['MASTER_ID'].to_list())]
        rhea_reactions['2_after_filtering_out_the_polymers'] = rhea_reactions['MASTER_ID'].isin(df_rhea_chebi['MASTER_ID'])
        
        print('Unique Rhea CLASS MASTER IDs without polymers, residues', len(df_rhea_chebi['MASTER_ID'].unique()))
        print('Unique Rhea chebi ids from class reactions without polymers, residues', len(df_rhea_chebi['chebi_id'].unique()))
        
        # Class Rhea to SwissLipids creates and processes the merged dataframe of Rhea reactions and Isomeric Subspecies
        # of ChEBIs in Rhea based on the SwissLipids ontology
        r2sl.get_df_swiss_lipids_chebi_rhea(df_swiss_lipids_chebi, df_rhea_chebi)

        # Step 1: Identify common columns for comparison
        common_cols = list(set(r2sl.df_swiss_lipids_chebi_rhea.columns) & set(df_rhea_chebi.columns))

        # Step 2: Find rows in df_rhea_chebi that are NOT in r2sl.df_swiss_lipids_chebi_rhea
        new_rows = df_rhea_chebi.merge(
            r2sl.df_swiss_lipids_chebi_rhea[common_cols],
            on=common_cols,
            how='left',
            indicator=True
        ).query('_merge == "left_only"').drop(columns=['_merge'])

        # Step 3: Append the new rows to r2sl.df_swiss_lipids_chebi_rhea
        r2sl.df_swiss_lipids_chebi_rhea = pd.concat([r2sl.df_swiss_lipids_chebi_rhea, new_rows], ignore_index=True)
        
        SLMs_in_rhea = r2sl.df_swiss_lipids_chebi_rhea['Lipid ID'].unique()
        print('Unique SLMs in Rhea:', len(SLMs_in_rhea))

        print('Total unique MASTER_ID:', len(r2sl.df_swiss_lipids_chebi_rhea['MASTER_ID'].unique()))
        print('Total unique chebi id:', len(r2sl.df_swiss_lipids_chebi_rhea['chebi_id'].unique()))
        print('Total unique rhea lipid id:', len(r2sl.df_swiss_lipids_chebi_rhea['Lipid ID'].unique()))

        print()
        master_ids_with_a_swisslipid = r2sl.df_swiss_lipids_chebi_rhea[r2sl.df_swiss_lipids_chebi_rhea['Lipid ID'].notna()]['MASTER_ID'].unique()
        df_master_ids_with_a_swisslipid = r2sl.df_swiss_lipids_chebi_rhea[r2sl.df_swiss_lipids_chebi_rhea['MASTER_ID'].isin(master_ids_with_a_swisslipid)]
        print('With SwissLipid unique MASTER_ID:', len(df_master_ids_with_a_swisslipid['MASTER_ID'].unique()))
        # print('With SwissLipid unique chebi id:', len(df_master_ids_with_a_swisslipid['chebi_id'].unique()))
        print('With SwissLipid unique rhea lipid id:', len(df_master_ids_with_a_swisslipid['Lipid ID'].unique()))

        print()
        master_ids_without_a_swisslipid = r2sl.df_swiss_lipids_chebi_rhea[~r2sl.df_swiss_lipids_chebi_rhea['MASTER_ID'].isin(master_ids_with_a_swisslipid)]['MASTER_ID'].unique()
        print('master_ids_without_a_swisslipid', len(master_ids_without_a_swisslipid))
        print('Example master_ids_without_a_swisslipid', random.sample(master_ids_without_a_swisslipid.tolist(), 3))

        master_id_star_compounds_star_is_not_lipid = df_master_ids_with_a_swisslipid[
            (
                (df_master_ids_with_a_swisslipid['smiles'].str.contains(r'\*', regex=True)) &
                (df_master_ids_with_a_swisslipid['Lipid ID'].isna())
            )
        ]['MASTER_ID'].unique()

        print('master_id_star_compounds_star_is_not_lipid', len(master_id_star_compounds_star_is_not_lipid))
        print('Examples master_id_star_compounds_star_is_not_lipid:', random.sample(master_id_star_compounds_star_is_not_lipid.tolist(), 3))

        # master_id_star_compound_is_lipid = df_master_ids_with_a_swisslipid[
        #     (
        #         (df_master_ids_with_a_swisslipid['smiles'].str.contains(r'\*', regex=True)) &
        #         (df_master_ids_with_a_swisslipid['rhea_lipid_id'].notna() | (df_master_ids_with_a_swisslipid['rhea_lipid_id'] == '') | (df_master_ids_with_a_swisslipid['rhea_lipid_id'] == 'NA'))
        #     )
        # ]['MASTER_ID'].unique()

        #df_master_id_star_compounds_star_is_lipid = df_master_ids_with_a_swisslipid[df_master_ids_with_a_swisslipid['MASTER_ID'].isin(master_id_star_compound_is_lipid)]
        df_master_id_star_compounds_star_is_lipid = df_master_ids_with_a_swisslipid[~df_master_ids_with_a_swisslipid['MASTER_ID'].isin(master_id_star_compounds_star_is_not_lipid)]
        print()
        print('df_master_id_star_compounds_star_is_lipid unique MASTER_ID:', len(df_master_id_star_compounds_star_is_lipid['MASTER_ID'].unique()))
        print('df_master_id_star_compounds_star_is_lipid unique chebi id:', len(df_master_id_star_compounds_star_is_lipid['chebi_id'].unique()))
        print('df_master_id_star_compounds_star_is_lipid unique rhea lipid id:', len(df_master_id_star_compounds_star_is_lipid['Lipid ID'].unique()))

        # Analyse the directed graph of SwissLipid ontology and get all isomeric subspecies per SLM in Rhea
        rhea_lipid_to_descendant_df = sl.get_lipid_to_descendant_df(SLMs_in_rhea)
        rhea_lipid_to_descendant_df.fillna('NA', inplace=True)
        # next df uses pos_descr_to_FA_list ->
        if no_curated_list_restrictions == False:
            print('CURATED')
            class_lipid_to_descendants_df, _ = sl.filter_curated_biologically_relevant_isomeric_subspecies_only(curated_fa_list_run=curated_fa_list_run)
            summary_lines.append(f'Total biologically human-relevant descendants of the class lipids identified in SwissLipids\t{len(class_lipid_to_descendants_df)}')
        elif no_curated_list_restrictions == True:
            class_lipid_to_descendants_df = rhea_lipid_to_descendant_df
        
        rhea_lipid_to_descendant_df_temp = rhea_lipid_to_descendant_df[rhea_lipid_to_descendant_df['Lipid ID'].isin(class_lipid_to_descendants_df['Lipid ID'])]

        r2sl.get_df_rhea_descendant(rhea_lipid_to_descendant_df_temp, sl.swisslipids)

        r2sl.df_rhea_descendant = r2sl.df_rhea_descendant[r2sl.df_rhea_descendant['MASTER_ID'].isin(df_master_id_star_compounds_star_is_lipid['MASTER_ID'])]
        print('Unique class lipid ids in Rhea:', len(r2sl.df_rhea_descendant['rhea_lipid_id'].unique()))
        print('Unique lipid isomeric subspecies descendants in Rhea:', len(r2sl.df_rhea_descendant['isomeric_subspecies_descendant_lipid_id'].unique()))

        df_star_compounds_isomeric_subspecies_na = r2sl.df_rhea_descendant[
            (
                (r2sl.df_rhea_descendant['smiles'].str.contains(r'\*', regex=True)) &
                (r2sl.df_rhea_descendant['isomeric_subspecies_descendant_lipid_id'].isna() | (r2sl.df_rhea_descendant['isomeric_subspecies_descendant_lipid_id'] == '') | (r2sl.df_rhea_descendant['isomeric_subspecies_descendant_lipid_id'] == 'NA'))
            )
        ]
        r2sl.df_rhea_descendant['isona'] = r2sl.df_rhea_descendant['isomeric_subspecies_descendant_lipid_id'].apply(lambda x: str(x).startswith('SLM'))
        
        print()
        print('df_star_compounds_isomeric_subspecies_na unique MASTER_ID:', len(df_star_compounds_isomeric_subspecies_na['MASTER_ID'].unique()))
        print('Example df_star_compounds_isomeric_subspecies_na', random.sample(df_star_compounds_isomeric_subspecies_na['MASTER_ID'].unique().tolist(), 3))
        print('df_star_compounds_isomeric_subspecies_na unique chebi id:', len(df_star_compounds_isomeric_subspecies_na['chebi_id'].unique()))
        print('df_star_compounds_isomeric_subspecies_na unique isomeric_subspecies_descendant_lipid_id:', len(df_star_compounds_isomeric_subspecies_na['isomeric_subspecies_descendant_lipid_id'].unique()))

        df_compounds_enumerated  = r2sl.df_rhea_descendant[~(r2sl.df_rhea_descendant['MASTER_ID'].isin(df_star_compounds_isomeric_subspecies_na['MASTER_ID'].unique()))]
        
        print()
        print('df_compounds_enumerated unique MASTER_ID:', len(df_compounds_enumerated['MASTER_ID'].unique()))
        print('df_compounds_enumerated unique chebi id:', len(df_compounds_enumerated['chebi_id'].unique()))


        df_star_compounds_enumerated = df_compounds_enumerated[
            (
                ~(df_compounds_enumerated['isomeric_subspecies_descendant_lipid_id'].isna() | (df_compounds_enumerated['isomeric_subspecies_descendant_lipid_id'] == '') | (df_compounds_enumerated['isomeric_subspecies_descendant_lipid_id'] == 'NA'))
            )
        ]
        print()
        print('df_star_compounds_enumerated unique MASTER_ID:', len( df_star_compounds_enumerated['MASTER_ID'].unique()))
        print('Example star_compounds_enumerated', random.sample(df_star_compounds_enumerated['MASTER_ID'].unique().tolist(), 3))
        print('df_star_compounds_enumerated unique chebi id:', len( df_star_compounds_enumerated['chebi_id'].unique()))
        print('df_star_compounds_enumerated unique rhea lipid id:', len( df_star_compounds_enumerated['rhea_lipid_id'].unique()))
        print('df_star_compounds_enumerated unique isomeric_subspecies_descendant_lipid_id:', len( df_star_compounds_enumerated['isomeric_subspecies_descendant_lipid_id'].unique()))

        print()

        ## Results Overview
        
        results_overview.dropna(subset=['Web-RInChIKey'], inplace=True)

        # Unique reactions
        print('--- Unique reactions ---')
        print('Enumerated reactions (with structure defined):', len(results_overview['Web-RInChIKey'].unique()))
        
        # MASTER_ID comparison
        master_ids_enumerated = set(results_overview['MASTER_ID'].to_list())
        print('Rhea master IDs used in enumeration', len(master_ids_enumerated))
        

        def extract_numeric_ids(df, col):
            """Extract all numeric IDs from a column, returned as a set of integers."""
            numbers = re.findall(r'\b\d+\b', ' '.join(df[col].dropna().astype(str).values))
            return set(map(int, numbers))

        def extract_slm_ids(df, col):
            """Extract all SLM IDs from a column."""
            return set(re.findall(r'SLM:\d+', ' '.join(df[col].dropna().astype(str).values)))

        # ChEBI ID comparison
        chebis = extract_numeric_ids(results_overview, 'chebi_equation')
        print("\n--- ChEBI ID ---")
        print(len(chebis))

        # SLM ID comparison
        slms = extract_slm_ids(results_overview, 'swisslipids_equation')
        print("\n--- SLM ID ---")
        print(len(slms))
        master_ids_enumerated = set([int(i) for i in master_ids_enumerated])
        compounds_enumerated = set(df_star_compounds_enumerated['MASTER_ID'].to_list())
        print('Only enumerated SL', len(set(df_star_compounds_enumerated['MASTER_ID'].to_list())-master_ids_enumerated))
        print(compounds_enumerated-master_ids_enumerated)
        print('Only in reactions', len(master_ids_enumerated-compounds_enumerated))
        print(master_ids_enumerated-compounds_enumerated)
        print('Both', len(set(results_overview['MASTER_ID'].to_list()).intersection(set(df_star_compounds_enumerated['MASTER_ID'].to_list()))))
        print()
        df_reaction_enumerated =  df_star_compounds_enumerated[df_star_compounds_enumerated['MASTER_ID'].isin(set(results_overview['MASTER_ID'].to_list()))]
        print('df_reaction_enumerated unique MASTER_ID:', len(df_reaction_enumerated['MASTER_ID'].unique()))
        print('df_reaction_enumerated unique chebi id:', len(df_reaction_enumerated['chebi_id'].unique()))
        print('df_reaction_enumerated unique rhea lipid id:', len(df_reaction_enumerated['rhea_lipid_id'].unique()))
        print('df_reaction_enumerated unique isomeric_subspecies_descendant_lipid_id:', len(df_reaction_enumerated['isomeric_subspecies_descendant_lipid_id'].unique()))

        df_input_for_graph = df_star_compounds_enumerated[df_star_compounds_enumerated['MASTER_ID'].isin(set(df_star_compounds_enumerated['MASTER_ID'].to_list())-set(results_overview['MASTER_ID'].to_list()))]
        #df_intput_for_graph.drop_duplicates(subset=['MASTER_ID', 'chebi_id', 'rhea_lipid_id'], inplace=True)
        df_input_for_graph.to_csv('input_for_graph.tsv', sep='\t', index=False)