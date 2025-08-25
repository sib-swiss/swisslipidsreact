import pandas as pd
import glob
import os
from datetime import datetime

from pyrheadb.RheaDB import RheaDB

from .SwissLipids import SwissLipids
from .RheaToSwisslipidsDf import RheaToSwisslipidsDf

def run_pipeline(curated_fa_list_run=True, output_dir=None, no_curated_list_restrictions=True, rheaid=None):

    # determine the base directory
    if output_dir is None:
        output_dir = os.getcwd()
    else:
        os.makedirs(output_dir, exist_ok=True)

    # generate a timestamped filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_results_file = f'{timestamp}_summary_results.tsv'

    # ---------- Load and Process RheaDB Data ----------

    # Top summary lines
    summary_lines = []

    # First read SwissLipids
    sl = SwissLipids(output_dir=output_dir, timestamp=timestamp)
    sl.read_swisslipids_from_file()

    # Map SwissLipids to ChEBI, and merge with Rhea data
    df_swiss_lipids_chebi = sl.get_only_chebi_sl_df()

    summary_lines.append(f'# swiss lipids * chebi\t{len(df_swiss_lipids_chebi)}')

    rdb = RheaDB()
    rheadf = rdb.rhea_reaction_long_format_smiles_chebi
    if rheaid:
        rheadf = rheadf[rheadf['MASTER_ID']==rheaid]
    rhea_reactions = rdb.df_reactions

    rhea_non_residue_reactions = rhea_reactions[rhea_reactions['residue_rxn_flag']==False]
    rhea_class_reactions = rhea_reactions[rhea_reactions['Web-RInChIKey'].isna()]
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
    df_rhea_prolymer_reactions = rheadf[rheadf['class'] == 'POLYMER']
    df_rhea_chebi = rheadf[~rheadf['MASTER_ID'].isin(df_rhea_prolymer_reactions['MASTER_ID'].to_list())]
    rhea_reactions['2_after_filtering_out_the_polymers'] = rhea_reactions['MASTER_ID'].isin(df_rhea_chebi['MASTER_ID'])
    
    # Class Rhea to SwissLipids creates and processes the merged dataframe of Rhea reactions and Isomeric Subspecies
    # of ChEBIs in Rhea based on the SwissLipids ontology
    r2sl = RheaToSwisslipidsDf(output_dir=output_dir, timestamp=timestamp)
    r2sl.get_df_swiss_lipids_chebi_rhea(df_swiss_lipids_chebi, df_rhea_chebi)
    SLMs_in_rhea = r2sl.df_swiss_lipids_chebi_rhea['Lipid ID'].unique()

    MASTER_IDs = r2sl.df_swiss_lipids_chebi_rhea['MASTER_ID'].unique()
    summary_lines.append(f'Class reaction Rhea IDs without polymers (lipid + non-lipid)\t{len(MASTER_IDs)}')

    rhea_reactions['3_keeping_only_reactions_with_swisslipids'] = rhea_reactions['MASTER_ID'].isin(MASTER_IDs)

    rhea_reactions['rxnsmiles_I'] = rhea_reactions['rxnsmiles'].apply(lambda x: x.replace('[1*]','C').replace('At','C').replace('[2*]','C').replace('*','C'))

    # def get_attention_guided_atom_map_error_wrap(mapped_rxn):
    #     try:
    #         return r2sl.rxn_mapper.get_attention_guided_atom_maps([mapped_rxn])[0]['mapped_rxn']
    #     except:
    #         return '[C:1]>>[C:1]'
    # rhea_reactions['rxnsmiles_I_mapped'] = rhea_reactions['rxnsmiles_I'].progress_apply(lambda x: get_attention_guided_atom_map_error_wrap(x))
    # rhea_reactions['bond_changes'] = rhea_reactions['rxnsmiles_I_mapped'].progress_apply(lambda x: len(r2sl.mapped_reaction_to_report(x)['bond_changes']))
    # rhea_reactions.to_csv('rhea_reactions_bond_changes.tsv', sep='\t', index=False)
    # exit()

    df_bc = pd.read_csv('rhea_reactions_bond_changes.tsv', sep='\t')
    bond_changes_lookup = dict(zip(df_bc['MASTER_ID'], df_bc['bond_changes']))

    # Analyse the directed graph of SwissLipid ontology and get all isomeric subspecies per SLM in Rhea
    rhea_lipid_to_descendant_df = sl.get_lipid_to_descendant_df(SLMs_in_rhea)
    summary_lines.append(f'# swiss lipids * chebi * rhea\t{len(r2sl.df_swiss_lipids_chebi_rhea)}')
    summary_lines.append(f'isomeric subspecies descendants of lipids * Rhea\t{len(rhea_lipid_to_descendant_df)}')
    summary_lines.append(f'unique Lipid ID isomeric subspecies descendants of lipids in Rhea\t{len(set(rhea_lipid_to_descendant_df["isomeric_subspecies_descendant_lipid_id"]))}')

    # next df uses pos_descr_to_FA_list ->
    if no_curated_list_restrictions == False:
        class_lipid_to_descendants_df, all_lipids_considered = sl.filter_curated_biologically_relevant_isomeric_subspecies_only(curated_fa_list_run=curated_fa_list_run)
        summary_lines.append(f'Total biologically human-relevant descendants of the class lipids identified in SwissLipids\t{len(class_lipid_to_descendants_df)}')
        rhea_lipid_to_descendant_df_not_filtered=rhea_lipid_to_descendant_df[~rhea_lipid_to_descendant_df['Lipid ID'].isin(all_lipids_considered)]
        class_lipid_to_descendants_df = pd.concat([class_lipid_to_descendants_df, rhea_lipid_to_descendant_df_not_filtered]) # adding those lipids that were not in the list with positions
    elif no_curated_list_restrictions == True:
        class_lipid_to_descendants_df = rhea_lipid_to_descendant_df
    
    rhea_lipid_to_descendant_df_temp = rhea_lipid_to_descendant_df[rhea_lipid_to_descendant_df['Lipid ID'].isin(class_lipid_to_descendants_df['Lipid ID'])]

    Rhea_x_swisslipid_isomeric_subspecies, \
        Unique_swisslipid_isomeric_subspecies_in_Rhea, \
            = r2sl.get_df_rhea_descendant(rhea_lipid_to_descendant_df_temp, sl.swisslipids)

    rhea_reactions['4_after_getting_df_descendants'] = rhea_reactions['MASTER_ID'].isin(r2sl.df_rhea_descendant['MASTER_ID'])

    # Function save_results generated the combinations and saves the results
    total_gen_attempted, \
    MASTER_ID_for_specific_fatty_acids_before_filtering_out_the_unbalanced, \
    num_reactions_to_check_for_balance, \
    MASTER_ID_for_specific_fatty_acids_after_filtering_out_the_unbalanced, \
    num_reactions_after_filtering_out_the_unbalanced = r2sl.save_results(rheadf, r2sl.df_rhea_descendant, rhea_reactions, f'{timestamp}_enumerated_reactions.tsv', bond_changes_lookup)

    stats_dict = {
        'Rhea_x_swisslipid_isomeric_subspecies':  Rhea_x_swisslipid_isomeric_subspecies, 
        'Unique_swisslipid_isomeric_subspecies_in_Rhea': Unique_swisslipid_isomeric_subspecies_in_Rhea,
        'Total MASTER_IDs with at least one SwissLipid descendant': total_gen_attempted, 
        'Total MASTER_IDs with SwissLipid in reactants and products of reaction':MASTER_ID_for_specific_fatty_acids_before_filtering_out_the_unbalanced,
        'num_reactions_to_check_for_balance': num_reactions_to_check_for_balance,
        'num_reactions_after_filtering_out_the_unbalanced':  num_reactions_after_filtering_out_the_unbalanced,
        'Master IDs after balance': MASTER_ID_for_specific_fatty_acids_after_filtering_out_the_unbalanced
    }

    # Save all to .tsv
    with open(os.path.join(output_dir, summary_results_file), 'w') as f:
        for line in summary_lines:
            f.write(f'{line}\n')
        for key, value in stats_dict.items():
            f.write(f'{key}\t{value}\n')

    print(f'Summary saved to {summary_results_file}')

