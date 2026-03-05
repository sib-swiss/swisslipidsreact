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

        # Read SwissLipids files.
        sl = SwissLipids(output_dir=output_dir)
        sl.read_swisslipids_from_file()

        print('Unique SLMs', len(sl.swisslipids['Lipid ID'].unique()))
        print('Unique SLMs isomeric subspecies', len(sl.swisslipids[sl.swisslipids['Level'] == 'Isomeric subspecies']['Lipid ID'].unique()))

        r2sl = RheaToSwisslipidsDf()

        # df_slm2chebi is a map of SLM ID to CHEBI ID:
        # Lipid ID | CHEBI | chebi_id
        # NB: 'CHEBI' is original value from lipids.tsv, 'chebi_id' has expanded pipe-separated 'CHEBI' values.
        df_slm2chebi = sl.build_df_slm2chebi()

        summary_lines.append(f'# swiss lipids * chebi\t{len(df_slm2chebi)}')
        print('Unique ChEBI IDs in SwissLipids:', len(df_slm2chebi['CHEBI'].unique()))
        print('Unique SL IDs with ChEBI:', len(df_slm2chebi['Lipid ID'].unique()))

        # Read Rhea files.
        rdb = RheaDB()

        # df_rhea2participants has one line per Rhea participant:
        # MASTER_ID | reaction_side | chebiid | smiles | inchi | inchikey | inchikey14L | stoich_coef
        # NB: 'chebiid' also contains POLYMER:x, but for GENERICs it contains the residues' CHEBI IDs.
        df_rhea2participants = rdb.rhea_reaction_long_format_smiles_chebi

        # df_rhea has one line per Rhea forward reaction:
        # MASTER_ID | reaction_participant_names | rheaid | rxnsmiles_no_residue_correction | chebi_equation | residue_rxn_flag | rxnsmiles | rxn_smiles_halogen | rxn_smiles_no_A_AH | class_reaction_flag | polymer_reaction_flag | RInChI | Web-RInChIKey
        df_rhea = rdb.df_reactions

        # Get the subset of Rhea reactions with class compounds (compounds with R-groups).
        # To do this we must first remove reactions with Rhea generic compounds,
        # because like those with class compounds they have no Web-RInChIKey.
        df_rhea_without_generics  = df_rhea[df_rhea['residue_rxn_flag'] == False]
        df_rhea_without_rinchikey = df_rhea[df_rhea['Web-RInChIKey'].isna()]
        df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'].isin(df_rhea_without_generics['MASTER_ID'])]
        df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'].isin(df_rhea_without_rinchikey['MASTER_ID'])]
        df_rhea = df_rhea.copy()
        df_rhea['1_df_rhea'] = df_rhea['MASTER_ID'].isin(df_rhea2participants['MASTER_ID'])

        #### SWISS LIPIDS CORRECTIONS ####
        # In swisslipids ChEBI:58342 corresponds to Fatty acyl-CoAs, however, in
        # ChEBI it is acyl-CoA, and in swisslipids CHEBI:77636 does not exist
        # Yes, it is confusing, thing is fatty acyl is a child of acyl,
        # but acyl is not always a lipid.
        # However acyl was used since the beginning of times to describe the fatty acyl reactions.
        # ->
        # I will change all CHEBI:77636 to CHEBI:58342 in the input to address this
        df_rhea2participants.loc[df_rhea2participants['chebiid'] == 'CHEBI:77636', 'chebiid'] = 'CHEBI:58342'

        # Replacing 2-monolysocardiolipin(2−) and 1-monolysocardiolipin(2−) with monolysocardiolipin(2−), since it is the enumerated version in SwissLipids
        df_rhea2participants.loc[df_rhea2participants['chebiid'] == 'CHEBI:64743', 'chebiid'] = 'CHEBI:167057'
        df_rhea2participants.loc[df_rhea2participants['chebiid'] == 'CHEBI:65092', 'chebiid'] = 'CHEBI:167057'

        # Replacing 2,3-diacyl-sn-glycerol with 1,2-diglyceride
        df_rhea2participants.loc[df_rhea2participants['chebiid'] == 'CHEBI:75524', 'chebiid'] = 'CHEBI:49172'

        # changing carboxilate to fatty acid for RHEA:34359 and RHEA:36195
        df_rhea2participants.loc[(df_rhea2participants['chebiid'] == 'CHEBI:29067') & (df_rhea2participants['MASTER_ID'] == 34359), 'chebiid'] = 'CHEBI:28868'
        df_rhea2participants.loc[(df_rhea2participants['chebiid'] == 'CHEBI:29067') & (df_rhea2participants['MASTER_ID'] == 36195), 'chebiid'] = 'CHEBI:28868'

        # changing aliphatic alcool to fatty alcohol for RHEA:59388
        df_rhea2participants.loc[(df_rhea2participants['chebiid'] == 'CHEBI:2571') & (df_rhea2participants['MASTER_ID'] == 59388), 'chebiid'] = 'CHEBI:24026'

        # changing a long chain fatty alcohol to fatty alcohol for RHEA:38443 since info about the chain length will come from the wax ester
        df_rhea2participants.loc[(df_rhea2participants['chebiid'] == 'CHEBI:17135') & (df_rhea2participants['MASTER_ID'] == 38443), 'chebiid'] = 'CHEBI:24026'

        # changing sterol to cholesterol since we only put cholesterol esters into the enumeration strategy
        # have to hardcode the SMILES for cholesterol (copied from SwissLipids)
        df_rhea2participants.loc[df_rhea2participants['chebiid'] == 'CHEBI:15889', 'smiles'] = 'CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C'
        df_rhea2participants.loc[df_rhea2participants['chebiid'] == 'CHEBI:15889', 'chebiid'] = 'CHEBI:16113'

        # Helper function to split a Rhea compound accession into prefix and numeric ID.
        def __split_compound_accession(compound):
            """
            Splits a Rhea compound accession (e.g., "CHEBI:1234" or "POLYMER:1234)
            into prefix and ID. Returns (prefix, ID) or (None, None) on failure.
            """
            try:
                prefix, number = compound.split(':')
                return pd.Series([prefix, int(number)])
            except Exception:
                return pd.Series([None, None])
        
        # Here we add a column 'id_prefix' to filter out the 'POLYMER:'.    
        df_rhea2participants[['id_prefix', 'chebi_id']] = df_rhea2participants['chebiid'].apply(__split_compound_accession)
        df_rhea2participants.dropna(subset=['chebi_id'], inplace=True)
        df_rhea2participants.loc[:, 'chebi_id'] = df_rhea2participants['chebi_id'].astype(int)

        # Remove reactions with polymers.
        l_rhea_with_polymers = df_rhea2participants.loc[df_rhea2participants['id_prefix'] == 'POLYMER', 'MASTER_ID'].tolist()
        df_rhea2participants = df_rhea2participants[~df_rhea2participants['MASTER_ID'].isin(l_rhea_with_polymers)]
        df_rhea['2_after_filtering_out_the_polymers'] = df_rhea['MASTER_ID'].isin(df_rhea2participants['MASTER_ID'])
        
        print('Unique Rhea CLASS MASTER IDs without polymers, residues', len(df_rhea2participants['MASTER_ID'].unique()))
        print('Unique Rhea chebi ids from class reactions without polymers, residues', len(df_rhea2participants['chebi_id'].unique()))
        
        # Class RheaToSwissLipids merges Rhea and SwissLipids data.
        # df_rhea2participants2slm_class has one line per Rhea participant:
        # MASTER_ID | reaction_side | chebiid | smiles | inchi | inchikey | inchikey14L | stoich_coef | id_prefix | chebi_id | Lipid ID | CHEBI
        # NB: rhea.chebi_id = sl.chebi_id -> Lipid ID = SLM ID of rhea.chebi_id
        r2sl.build_df_rhea2participants2slm_class(df_rhea2participants, df_slm2chebi)

        # Step 1: Identify common columns for comparison
        common_cols = list(set(r2sl.df_rhea2participants2slm_class.columns) & set(df_rhea2participants.columns))

        # Step 2: Find rows in df_rhea2participants that are NOT in r2sl.df_rhea2participants2slm_class
        new_rows = df_rhea2participants.merge(
            r2sl.df_rhea2participants2slm_class[common_cols],
            on=common_cols,
            how='left',
            indicator=True
        ).query('_merge == "left_only"').drop(columns=['_merge'])

        # Step 3: Append the new rows to r2sl.df_rhea2participants2slm_class
        r2sl.df_rhea2participants2slm_class = pd.concat([r2sl.df_rhea2participants2slm_class, new_rows], ignore_index=True)
        
        SLMs_in_rhea = r2sl.df_rhea2participants2slm_class['Lipid ID'].unique()
        print('Unique SLMs in Rhea:', len(SLMs_in_rhea))

        print('Total unique MASTER_ID:',     len(r2sl.df_rhea2participants2slm_class['MASTER_ID'].unique()))
        print('Total unique chebi id:',      len(r2sl.df_rhea2participants2slm_class['chebi_id'].unique()))
        print('Total unique rhea lipid id:', len(r2sl.df_rhea2participants2slm_class['Lipid ID'].unique()))

        print()
        master_ids_with_a_swisslipid = r2sl.df_rhea2participants2slm_class[r2sl.df_rhea2participants2slm_class['Lipid ID'].notna()]['MASTER_ID'].unique()
        df_master_ids_with_a_swisslipid = r2sl.df_rhea2participants2slm_class[r2sl.df_rhea2participants2slm_class['MASTER_ID'].isin(master_ids_with_a_swisslipid)]
        print('With SwissLipid unique MASTER_ID:', len(df_master_ids_with_a_swisslipid['MASTER_ID'].unique()))
        # print('With SwissLipid unique chebi id:', len(df_master_ids_with_a_swisslipid['chebi_id'].unique()))
        print('With SwissLipid unique rhea lipid id:', len(df_master_ids_with_a_swisslipid['Lipid ID'].unique()))

        print()
        master_ids_without_a_swisslipid = r2sl.df_rhea2participants2slm_class[~r2sl.df_rhea2participants2slm_class['MASTER_ID'].isin(master_ids_with_a_swisslipid)]['MASTER_ID'].unique()
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
        #         (df_master_ids_with_a_swisslipid['class_slm_id'].notna() | (df_master_ids_with_a_swisslipid['rhea_slm_id'] == '') | (df_master_ids_with_a_swisslipid['rhea_slm_id'] == 'NA'))
        #     )
        # ]['MASTER_ID'].unique()

        #df_master_id_star_compounds_star_is_lipid = df_master_ids_with_a_swisslipid[df_master_ids_with_a_swisslipid['MASTER_ID'].isin(master_id_star_compound_is_lipid)]
        df_master_id_star_compounds_star_is_lipid = df_master_ids_with_a_swisslipid[~df_master_ids_with_a_swisslipid['MASTER_ID'].isin(master_id_star_compounds_star_is_not_lipid)]
        print()
        print('df_master_id_star_compounds_star_is_lipid unique MASTER_ID:', len(df_master_id_star_compounds_star_is_lipid['MASTER_ID'].unique()))
        print('df_master_id_star_compounds_star_is_lipid unique chebi id:', len(df_master_id_star_compounds_star_is_lipid['chebi_id'].unique()))
        print('df_master_id_star_compounds_star_is_lipid unique rhea lipid id:', len(df_master_id_star_compounds_star_is_lipid['Lipid ID'].unique()))

        # Analyse the directed graph of SwissLipid ontology and get all isomeric subspecies per SLM in Rhea
        df_rhea_slm_class2iso = sl.build_df_slm_class2iso(SLMs_in_rhea)
        df_rhea_slm_class2iso.fillna('NA', inplace=True)
        # next df uses pos_descr_to_FA_list ->
        if no_curated_list_restrictions == False:
            print('CURATED')
            df_rhea_slm_class2iso_filtered, _ = sl.filter_curated_biologically_relevant_isomeric_subspecies_only(curated_fa_list_run=curated_fa_list_run)
            summary_lines.append(f'Total biologically human-relevant descendants of the class lipids identified in SwissLipids\t{len(df_rhea_slm_class2iso_filtered)}')
        elif no_curated_list_restrictions == True:
            df_rhea_slm_class2iso_filtered = df_rhea_slm_class2iso
        
        df_rhea_slm_class2iso_temp = df_rhea_slm_class2iso[df_rhea_slm_class2iso['Lipid ID'].isin(df_rhea_slm_class2iso_filtered['Lipid ID'])]

        r2sl.build_df_rhea2participants2slm_class2iso(df_rhea_slm_class2iso_temp, sl.swisslipids)

        r2sl.df_rhea2participants2slm_class2iso = r2sl.df_rhea2participants2slm_class2iso[r2sl.df_rhea2participants2slm_class2iso['MASTER_ID'].isin(df_master_id_star_compounds_star_is_lipid['MASTER_ID'])]
        print('Unique class lipid ids in Rhea:', len(r2sl.df_rhea2participants2slm_class2iso['rhea_slm_id'].unique()))
        print('Unique lipid isomeric subspecies descendants in Rhea:', len(r2sl.df_rhea2participants2slm_class2iso['iso_slm_id'].unique()))

        df_star_compounds_isomeric_subspecies_na = r2sl.df_rhea2participants2slm_class2iso[
            (
                (r2sl.df_rhea2participants2slm_class2iso['smiles'].str.contains(r'\*', regex=True)) &
                (r2sl.df_rhea2participants2slm_class2iso['iso_slm_id'].isna() | (r2sl.df_rhea2participants2slm_class2iso['iso_slm_id'] == '') | (r2sl.df_rhea2participants2slm_class2iso['iso_slm_id'] == 'NA'))
            )
        ]
        r2sl.df_rhea2participants2slm_class2iso['isona'] = r2sl.df_rhea2participants2slm_class2iso['iso_slm_id'].apply(lambda x: str(x).startswith('SLM'))
        
        print()
        print('df_star_compounds_isomeric_subspecies_na unique MASTER_ID:', len(df_star_compounds_isomeric_subspecies_na['MASTER_ID'].unique()))
        print('Example df_star_compounds_isomeric_subspecies_na', random.sample(df_star_compounds_isomeric_subspecies_na['MASTER_ID'].unique().tolist(), 3))
        print('df_star_compounds_isomeric_subspecies_na unique chebi id:', len(df_star_compounds_isomeric_subspecies_na['chebi_id'].unique()))
        print('df_star_compounds_isomeric_subspecies_na unique iso_slm_id:', len(df_star_compounds_isomeric_subspecies_na['iso_slm_id'].unique()))

        df_compounds_enumerated  = r2sl.df_rhea2participants2slm_class2iso[~(r2sl.df_rhea2participants2slm_class2iso['MASTER_ID'].isin(df_star_compounds_isomeric_subspecies_na['MASTER_ID'].unique()))]
        
        print()
        print('df_compounds_enumerated unique MASTER_ID:', len(df_compounds_enumerated['MASTER_ID'].unique()))
        print('df_compounds_enumerated unique chebi id:', len(df_compounds_enumerated['chebi_id'].unique()))


        df_star_compounds_enumerated = df_compounds_enumerated[
            (
                ~(df_compounds_enumerated['iso_slm_id'].isna() | (df_compounds_enumerated['iso_slm_id'] == '') | (df_compounds_enumerated['iso_slm_id'] == 'NA'))
            )
        ]
        print()
        print('df_star_compounds_enumerated unique MASTER_ID:', len( df_star_compounds_enumerated['MASTER_ID'].unique()))
        print('Example star_compounds_enumerated', random.sample(df_star_compounds_enumerated['MASTER_ID'].unique().tolist(), 3))
        print('df_star_compounds_enumerated unique chebi id:', len( df_star_compounds_enumerated['chebi_id'].unique()))
        print('df_star_compounds_enumerated unique rhea lipid id:', len( df_star_compounds_enumerated['rhea_slm_id'].unique()))
        print('df_star_compounds_enumerated unique iso_slm_id:', len( df_star_compounds_enumerated['iso_slm_id'].unique()))

        print()

        ## Results Overview
        
        results_overview.dropna(subset=['Web-RInChIKey'], inplace=True)

        # Unique reactions
        print('--- Unique reactions ---')
        print('Enumerated reactions (with structure defined):', len(results_overview['Web-RInChIKey'].unique()))
        
        # MASTER_ID comparison
        master_ids_enumerated = set(results_overview['MASTER_ID'].tolist())
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
        compounds_enumerated = set(df_star_compounds_enumerated['MASTER_ID'].tolist())
        print('Only enumerated SL', len(set(df_star_compounds_enumerated['MASTER_ID'].tolist())-master_ids_enumerated))
        print(compounds_enumerated-master_ids_enumerated)
        print('Only in reactions', len(master_ids_enumerated-compounds_enumerated))
        print(master_ids_enumerated-compounds_enumerated)
        print('Both', len(set(results_overview['MASTER_ID'].tolist()).intersection(set(df_star_compounds_enumerated['MASTER_ID'].tolist()))))
        print()
        df_reaction_enumerated =  df_star_compounds_enumerated[df_star_compounds_enumerated['MASTER_ID'].isin(set(results_overview['MASTER_ID'].tolist()))]
        print('df_reaction_enumerated unique MASTER_ID:', len(df_reaction_enumerated['MASTER_ID'].unique()))
        print('df_reaction_enumerated unique chebi IDs:', len(df_reaction_enumerated['chebi_id'].unique()))
        print('df_reaction_enumerated unique rhea SLM IDs:', len(df_reaction_enumerated['rhea_slm_id'].unique()))
        print('df_reaction_enumerated unique isomeric SLM IDs:', len(df_reaction_enumerated['iso_slm_id'].unique()))

        df_input_for_graph = df_star_compounds_enumerated[df_star_compounds_enumerated['MASTER_ID'].isin(set(df_star_compounds_enumerated['MASTER_ID'].tolist())-set(results_overview['MASTER_ID'].tolist()))]
        #df_intput_for_graph.drop_duplicates(subset=['MASTER_ID', 'chebi_id', 'rhea_slm_id'], inplace=True)
        df_input_for_graph.to_csv('input_for_graph.tsv', sep='\t', index=False)
