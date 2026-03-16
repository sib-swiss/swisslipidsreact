import logging
import pandas as pd
import glob
import os
import re

from datetime import datetime

from pyrheadb.RheaDB import RheaDB

from .utils import DEBUG, debug_df_first_row
from .SwissLipids import SwissLipids
from .RheaToSwisslipidsDf import RheaToSwisslipidsDf

logger = logging.getLogger(__name__)

def run_pipeline(output_dir=None, filter_fa="curated", filter_rhea=False, rhea_id=None, rhea_version=None, test=False):

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

    # Read SwissLipids files.
    sl = SwissLipids(output_dir=output_dir, timestamp=timestamp)
    sl.read_swisslipids_from_file(filter_fa, test)

    # df_slm2chebi maps an SLM ID to a CHEBI ID:
    # Lipid ID | CHEBI | chebi_id
    # NB: 'CHEBI' is original value from lipids.tsv, 'chebi_id' has expanded pipe-separated 'CHEBI' values.
    df_slm2chebi = sl.build_df_slm2chebi()
    logger.info( "LENGTH: %8d df_slm2chebi", len(df_slm2chebi) )
    logger.debug( "FORMAT: df_slm2chebi\n%s", debug_df_first_row(df_slm2chebi) )
    if DEBUG > 1:
        df_slm2chebi.to_csv(os.path.join(output_dir, 'DEBUG_df_slm2chebi.tsv'), sep="\t", header=True, index=False)
    summary_lines.append(f'df_slm2chebi\t{len(df_slm2chebi)}')

    if filter_rhea == True:
        # Get the unique SLM lipid classes of all isomeric subspecies.
        df_iso_class = sl.get_isomeric_subspecies_class()
        logger.info( "LENGTH: %8d df_iso_class", len(df_iso_class) )
        logger.debug( "FORMAT: df_iso_class\n%s", debug_df_first_row(df_iso_class) )
        if DEBUG > 1:
            df_iso_class.to_csv(os.path.join(output_dir, 'DEBUG_df_iso_class.tsv'), sep="\t", header=True, index=False)

        # Get the ChEBI ID for each SLM ID.
        df_iso_class2chebi = df_iso_class.merge(df_slm2chebi, left_on='Lipid class*', right_on='Lipid ID', how='inner')
        logger.info( "LENGTH: %8d df_iso_class2chebi", len(df_iso_class2chebi) )
        logger.debug( "FORMAT: df_iso_class2chebi\n%s", debug_df_first_row(df_iso_class2chebi) )
        if DEBUG > 1:
            df_iso_class2chebi.to_csv(os.path.join(output_dir, 'DEBUG_df_iso_class2chebi.tsv'), sep="\t", header=True, index=False)

        # Get the list of ChEBIs.
        l_iso_class_chebi = df_iso_class2chebi['chebi_id'].tolist()
        logger.info( "Statistics: %8d l_iso_class_chebi", len(l_iso_class_chebi ) )
    
    # Read the Rhea FTP files.
    rdb = RheaDB(rhea_version=rhea_version)
    rhea_version = rdb.rhea_db_version
    logger.info( "Loaded files of Rhea release %s", rhea_version )

    # df_rhea has one line per Rhea forward reaction:
    # MASTER_ID | reaction_participant_names | rheaid | rxnsmiles_no_residue_correction | chebi_equation | residue_rxn_flag | rxnsmiles | rxn_smiles_halogen | rxn_smiles_no_A_AH | class_reaction_flag | polymer_reaction_flag | RInChI | Web-RInChIKey
    df_rhea = rdb.df_reactions
    logger.info( "LENGTH: %8d df_rhea", len(df_rhea) )
    logger.debug( "FORMAT: df_rhea\n%s", debug_df_first_row(df_rhea) )
    if DEBUG > 1:
        df_rhea.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea.tsv'), sep="\t", header=True, index=False)
    logger.info( "Statistics: %8d Rhea IDs (unfiltered)", len(df_rhea) )
    
    # df_rhea2participants has one line per Rhea participant:
    # MASTER_ID | reaction_side | chebiid | smiles | inchi | inchikey | inchikey14L | stoich_coef
    # NB: 'chebiid' also contains POLYMER:x, but for GENERICs it contains the residues' CHEBI IDs.
    df_rhea2participants = rdb.rhea_reaction_long_format_smiles_chebi
    logger.info( "LENGTH: %8d df_rhea2participants (unfiltered)", len(df_rhea2participants) )
    logger.debug( "FORMAT: df_rhea2participants\n%s", debug_df_first_row(df_rhea2participants) )
    if DEBUG > 1:
        df_rhea2participants.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea2participants.tsv'), sep="\t", header=True, index=False)
    if rhea_id:
        df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'] == rhea_id]
    
    # Get the subset of Rhea reactions with class compounds (compounds with R-groups).
    # To do this we must first remove reactions with Rhea generic compounds,
    # because like those with class compounds they have no Web-RInChIKey.
    df_rhea_without_generics  = df_rhea[df_rhea['residue_rxn_flag'] == False]
    df_rhea_without_rinchikey = df_rhea[df_rhea['Web-RInChIKey'].isna()]
    df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'].isin(df_rhea_without_generics['MASTER_ID'])]
    logger.info( "LENGTH: %8d df_rhea2participants (after filtering by generics)", len(df_rhea2participants) )
    logger.info( "Statistics: %8d Rhea IDs (after filtering by generics)", len(df_rhea2participants['MASTER_ID'].unique()) )
    df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'].isin(df_rhea_without_rinchikey['MASTER_ID'])]
    logger.info( "LENGTH: %8d df_rhea2participants (after filtering by Web-RInChIKey)", len(df_rhea2participants) )
    logger.info( "Statistics: %8d Rhea IDs (after filtering by Web-RInChIKey)", len(df_rhea2participants['MASTER_ID'].unique()) )
    df_rhea = df_rhea.copy()
    df_rhea['1_after_filtering_generics_and_rinchikey'] = df_rhea['MASTER_ID'].isin(df_rhea2participants['MASTER_ID'])
    
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
    logger.info( "LENGTH: %8d df_rhea2participants (after filtering by polymers)", len(df_rhea2participants) )
    logger.info( "Statistics: %8d Rhea IDs (after filtering by polymers)", len(df_rhea2participants['MASTER_ID'].unique()) )
    df_rhea['2.1_after_filtering_out_the_polymers'] = df_rhea['MASTER_ID'].isin(df_rhea2participants['MASTER_ID'])
    
    # Keep only reactions with at least one participant that was used as a class in the enumeration.
    # NB: We can only do this here, because it needs a clean column 'chebi_id'.
    if filter_rhea == True:
        l_rhea_used_in_enumeration = df_rhea2participants.loc[df_rhea2participants['chebi_id'].isin(l_iso_class_chebi), 'MASTER_ID'].unique().tolist()
        logger.info( "Statistics: %8d l_rhea_used_in_enumeration", len(l_rhea_used_in_enumeration) )

        df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'].isin(l_rhea_used_in_enumeration)]
        logger.info( "LENGTH: %8d df_rhea2participants (after filtering by SLM used in enumeration)", len(df_rhea2participants) )
        logger.info( "Statistics: %6d Potential reaction (after filtering by SLM used in enumeration)", len(df_rhea2participants['MASTER_ID'].unique()) )
        df_rhea['2.2_after_filtering_by_enumeration_classes'] = df_rhea['MASTER_ID'].isin(l_rhea_used_in_enumeration)
        
    # Class RheaToSwissLipids merges Rhea and SwissLipids data.
    r2sl = RheaToSwisslipidsDf(output_dir=output_dir, timestamp=timestamp)
    # df_rhea2participants2slm_class maps Rhea participants chebi_id to the corresponding SLM ID ('Lipid ID'):
    # MASTER_ID | reaction_side | chebiid | smiles | inchi | inchikey | inchikey14L | stoich_coef | id_prefix | chebi_id | Lipid ID | CHEBI
    r2sl.build_df_rhea2participants2slm_class(df_rhea2participants, df_slm2chebi)
    logger.info( "LENGTH: %8d df_rhea2participants2slm_class", len(r2sl.df_rhea2participants2slm_class) )
    logger.debug( "FORMAT: df_rhea2participants2slm_class\n%s", debug_df_first_row(r2sl.df_rhea2participants2slm_class) )
    if DEBUG > 1:
        r2sl.df_rhea2participants2slm_class.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea2participants2slm_class.tsv'), sep="\t", header=True, index=False)

    # Get the list of unique Rhea class SLM IDs.     
    SLMs_in_rhea = r2sl.df_rhea2participants2slm_class['Lipid ID'].unique()
    logger.info( "Statistics: %8d SLMs_in_rhea", len(SLMs_in_rhea) )
    MASTER_IDs = r2sl.df_rhea2participants2slm_class['MASTER_ID'].unique()
    summary_lines.append(f'Class reaction Rhea IDs without polymers (lipid + non-lipid)\t{len(MASTER_IDs)}')

    df_rhea['3_keeping_only_reactions_with_swisslipids'] = df_rhea['MASTER_ID'].isin(MASTER_IDs)

   # Compute and cache the bond changes (we use them to filter the enumeration results).
    # This takes app. 20min to compute, therefore we cache the result file in the rhea release folder.
    def __get_attention_guided_atom_map_error_wrap(mapped_rxn):
        try:
            return r2sl.rxn_mapper.get_attention_guided_atom_maps([mapped_rxn])[0]['mapped_rxn']
        except:
            return '[C:1]>>[C:1]' # return dummy atom mapped equation that will return 0 bond changes

    f_cache_rhea_bond_changes = os.path.join(rdb.rhea_db_version_location, "tsv/rhea-bond-changes-cache.tsv")
    if not os.path.exists(f_cache_rhea_bond_changes):
        logger.info( "Creating %s", f_cache_rhea_bond_changes )
        logger.info( "progress_apply: Generating rhea_reactions_bond_changes.tsv - Part 1" )
        # TODO: What is this? It seems to be fixing some SMILES?
        df_rhea['rxnsmiles_I'] = df_rhea['rxnsmiles'].apply(lambda x: x.replace('[1*]','C').replace('At','C').replace('[2*]','C').replace('*','C'))
        df_rhea['rxnsmiles_I_mapped'] = df_rhea['rxnsmiles_I'].progress_apply(lambda x: __get_attention_guided_atom_map_error_wrap(x))
        logger.info( "progress_apply: Generating rhea_reactions_bond_changes.tsv - Part 2" )
        df_rhea['bond_changes'] = df_rhea['rxnsmiles_I_mapped'].progress_apply(lambda x: len(r2sl.mapped_reaction_to_report(x)['bond_changes']))
        df_rhea.to_csv(f_cache_rhea_bond_changes, sep='\t', index=False)
        
    logger.info( "Reading %s", f_cache_rhea_bond_changes )
    df_rhea_bond_changes = pd.read_csv(f_cache_rhea_bond_changes, sep='\t')
    logger.info( "LENGTH: %8d df_rhea_bond_changes", len(df_rhea_bond_changes) )
    dict_rhea2bond_changes = dict(zip(df_rhea_bond_changes['MASTER_ID'], df_rhea_bond_changes['bond_changes']))

    # df_rhea_slm_class2iso maps a Rhea class SLM ID to its isomeric subspecies SLM IDs:
    # class_slm_id | iso_slm_id | Lipid ID | Level | Name | Lipid class* | Components* | SMILES (pH7.3) | CHEBI | sn1' | sn2' | ...
    # NB: iso_slm_id = Lipid ID
    df_rhea_slm_class2iso = sl.build_df_slm_class2iso(SLMs_in_rhea)
    logger.info( "LENGTH: %8d df_rhea_slm_class2iso", len(df_rhea_slm_class2iso) )
    logger.debug( "FORMAT: df_rhea_slm_class2iso\n%s", debug_df_first_row(df_rhea_slm_class2iso) )
    if DEBUG > 1:
        df_rhea_slm_class2iso.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea_slm_class2iso.tsv'), sep="\t", header=True, index=False)
    summary_lines.append(f'r2sl.df_rhea2participants2slm_class\t{len(r2sl.df_rhea2participants2slm_class)}')
    summary_lines.append(f'df_rhea_slm_class2iso\t{len(df_rhea_slm_class2iso)}')
    summary_lines.append(f'df_rhea_slm_class2iso unique iso_slm_id\t{len(set(df_rhea_slm_class2iso["iso_slm_id"]))}')
    
    df_rhea_slm_class2iso_filtered = df_rhea_slm_class2iso
    # Filter fatty acids.
    if filter_fa == "curated":
        df_rhea_slm_class2iso_filtered, all_lipids_considered = sl.filter_fa_curated(test=test)
        summary_lines.append(f'df_rhea_slm_class2iso_filtered (biologically human-relevant)\t{len(df_rhea_slm_class2iso_filtered)}')
        df_rhea_slm_class2iso_not_filtered = df_rhea_slm_class2iso[~df_rhea_slm_class2iso['Lipid ID'].isin(all_lipids_considered)]
        logger.info( "LENGTH: %8d df_rhea_slm_class2iso_filtered (curated)", len(df_rhea_slm_class2iso_filtered) )
        logger.info( "LENGTH: %8d df_rhea_slm_class2iso_not_filtered (curated)", len(df_rhea_slm_class2iso_not_filtered) )
        if DEBUG > 1:
            df_rhea_slm_class2iso_not_filtered.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea_slm_class2iso_not_filtered.tsv'), sep="\t", header=True, index=False)
        # Add those lipids that were not in the list with positions.
        df_rhea_slm_class2iso_filtered = pd.concat([df_rhea_slm_class2iso_filtered, df_rhea_slm_class2iso_not_filtered])
        
    # FORMAT: df_rhea_slm_class2iso_filtered = df_rhea_slm_class2iso
    logger.info( "LENGTH: %8d df_rhea_slm_class2iso_filtered", len(df_rhea_slm_class2iso_filtered) )
    logger.debug( "FORMAT: df_rhea_slm_class2iso_filtered\n%s", debug_df_first_row(df_rhea_slm_class2iso_filtered) )

    df_rhea_slm_class2iso_temp = df_rhea_slm_class2iso[df_rhea_slm_class2iso['Lipid ID'].isin(df_rhea_slm_class2iso_filtered['Lipid ID'])]
    logger.info( "LENGTH: %8d df_rhea_slm_class2iso_temp", len(df_rhea_slm_class2iso_temp) )

    Rhea_x_swisslipid_isomeric_subspecies, \
        Unique_swisslipid_isomeric_subspecies_in_Rhea, \
            = r2sl.build_df_rhea2participants2slm_class2iso(df_rhea_slm_class2iso_temp, sl.swisslipids)

    logger.info( "Statistics: %8d Rhea IDs (after mapping class to iso)", len(r2sl.df_rhea2participants2slm_class2iso['MASTER_ID'].unique()) )
    df_rhea['4_after_mapping_class2iso'] = df_rhea['MASTER_ID'].isin(r2sl.df_rhea2participants2slm_class2iso['MASTER_ID'])

    # Function save_results generated the combinations and saves the results
    total_gen_attempted, \
    MASTER_ID_for_specific_fatty_acids_before_filtering_out_the_unbalanced, \
    num_reactions_to_check_for_balance, \
    MASTER_ID_for_specific_fatty_acids_after_filtering_out_the_unbalanced, \
    num_reactions_after_filtering_out_the_unbalanced = r2sl.save_results(df_rhea2participants, df_rhea, f'{timestamp}_enumerated_reactions.tsv', dict_rhea2bond_changes)

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

    logger.info( "Summary saved to %s\n", summary_results_file )

    # df_rhea['x_...'] contains a Pandas boolean Series. Use sum() to get number of True.
    pattern = r'^\d+_'
    for key, value in df_rhea.items():
        if re.match(pattern, key):
            print(value.sum(), " - ", key)
