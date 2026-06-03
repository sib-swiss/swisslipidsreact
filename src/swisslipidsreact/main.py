import logging
import pandas as pd
import glob
import os
import re

from datetime import datetime

from pyrheadb.RheaDB import RheaDB

from .utils import *
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
    
    # ---------- Load and Process RheaDB Data ----------

    # Read SwissLipids files.
    sl = SwissLipids(output_dir=output_dir, timestamp=timestamp)
    sl.read_swisslipids_from_file(filter_fa, test)
    # df_slm2chebi maps an SLM ID to a CHEBI ID:
    # slm_id | chebi_id
    df_slm2chebi = sl.build_df_slm2chebi()
    
    # Read the Rhea FTP files.
    rdb = RheaDB(rhea_version=rhea_version)
    rhea_version = rdb.rhea_db_version
    logger.info( "Loaded files of Rhea release %s", rhea_version )

    # df_rhea has one line per Rhea forward reaction:
    # MASTER_ID | reaction_participant_names | rheaid | rxnsmiles_no_residue_correction | chebi_equation | residue_rxn_flag | rxnsmiles | rxn_smiles_halogen | rxn_smiles_no_A_AH | class_reaction_flag | polymer_reaction_flag | RInChI | Web-RInChIKey
    df_rhea = rdb.df_reactions
    add_statistics(logger, "rhea", "%8d Rhea IDs (unfiltered)\n" % len(df_rhea))
    if DEBUG > 1:
        df_rhea.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea.tsv'), sep="\t", header=True, index=False)
    
    # df_rhea2participants has one line per Rhea participant:
    # MASTER_ID | reaction_side | chebiid | smiles | inchi | inchikey | inchikey14L | stoich_coef
    # NB: 'chebiid' also contains POLYMER:x, but for GENERICs it contains the residues' CHEBI IDs.
    df_rhea2participants = rdb.rhea_reaction_long_format_smiles_chebi
    add_statistics(logger, "participants", "%8d Rhea participants (unfiltered)\n" % len(df_rhea2participants))
    if DEBUG > 1:
        df_rhea2participants.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea2participants.tsv'), sep="\t", header=True, index=False)
    if rhea_id:
        df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'] == rhea_id]
        
    # Get the subset of Rhea reactions with class compounds (compounds with R-groups).
    # 1. Remove reactions with Web-RInChIKey.
    df_rhea_without_rinchikey = df_rhea[df_rhea['Web-RInChIKey'].isna()]
    df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'].isin(df_rhea_without_rinchikey['MASTER_ID'])]
    add_statistics(logger, "participants", "%8d Rhea participants (after filtering by Web-RInChIKey)\n" % len(df_rhea2participants))
    add_statistics(logger, "rhea", "%8d Rhea IDs (after filtering by Web-RInChIKey)\n" % len(df_rhea2participants['MASTER_ID'].unique()))
    # 2. Remove reactions with Rhea generic compounds
    df_rhea_without_generics  = df_rhea[df_rhea['residue_rxn_flag'] == False]
    df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'].isin(df_rhea_without_generics['MASTER_ID'])]
    add_statistics(logger, "participants", "%8d Rhea participants (after filtering by generics)\n" % len(df_rhea2participants))
    add_statistics(logger, "rhea", "%8d Rhea IDs (after filtering by generics)\n" % len(df_rhea2participants['MASTER_ID'].unique()))    
    
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
    add_statistics(logger, "participants", "%8d Rhea participants (after filtering by polymers)\n" % len(df_rhea2participants))
    add_statistics(logger, "rhea", "%8d Rhea IDs (after filtering by polymers)\n" % len(df_rhea2participants['MASTER_ID'].unique()))
    
    # Class RheaToSwissLipids merges Rhea and SwissLipids data.
    r2sl = RheaToSwisslipidsDf(output_dir=output_dir, timestamp=timestamp)
    # df_rhea2participants2slm maps Rhea participants chebi_id to the corresponding slm_id:
    # MASTER_ID | reaction_side | chebiid | smiles | inchi | inchikey | inchikey14L | stoich_coef | id_prefix | chebi_id | slm_id | CHEBI
    r2sl.build_df_rhea2participants2slm(df_rhea2participants, df_slm2chebi)
    add_statistics(logger, "participants", "%8d Rhea participants (after mapping participant ChEBI ID to SLM ID)\n" % len(r2sl.df_rhea2participants2slm))
    add_statistics(logger, "rhea", "%8d Rhea IDs (after mapping participant ChEBI ID to SLM ID)\n" % len(r2sl.df_rhea2participants2slm['MASTER_ID'].unique()))
    if DEBUG > 1:
        r2sl.df_rhea2participants2slm.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea2participants2slm.tsv'), sep="\t", header=True, index=False)

    # Get the list of unique SwissLipids SLM IDs in Rhea.
    SLMs_in_rhea = r2sl.df_rhea2participants2slm['slm_id'].unique()
    # Check which parent classes are not in Rhea: some are candidates for new reactions?
    df_iso_parent_in_rhea     = (sl.df_slm_parent2iso[sl.df_slm_parent2iso['parent_slm_id'].isin(SLMs_in_rhea)]['parent_slm_id'].drop_duplicates())
    df_iso_parent_not_in_rhea = (sl.df_slm_parent2iso[~sl.df_slm_parent2iso['parent_slm_id'].isin(SLMs_in_rhea)]['parent_slm_id'].drop_duplicates())
    add_statistics(logger, "lipids", "-------- SwissLipids parent class SLMs IDs\n")
    add_statistics(logger, "lipids", "%8d SwissLipids parent class SLMs IDs\n" % len(set(sl.df_slm_parent2iso['parent_slm_id'])))
    add_statistics(logger, "lipids", "%8d SwissLipids parent class SLMs IDs in Rhea\n" % len(df_iso_parent_in_rhea))
    add_statistics(logger, "lipids", "%8d SwissLipids parent class SLMs IDs not in Rhea\n" % len(df_iso_parent_not_in_rhea))
    if DEBUG > 1:
        df_iso_parent_in_rhea.to_csv(os.path.join(output_dir, 'DEBUG_df_iso_parent_in_rhea.tsv'), sep="\t", header=True, index=False) 
        df_iso_parent_not_in_rhea.to_csv(os.path.join(output_dir, 'DEBUG_df_iso_parent_not_in_rhea.tsv'), sep="\t", header=True, index=False)
    
    # Keep only reactions with participants that have been used as a parent class in the enumeration.
    # NB: We can only do this here, because it needs a clean column 'chebi_id'.
    # - Option: EACH side has a participant that is a parent class SLMs IDs of isomeric subspecies.
    if filter_rhea == "two-sides":
        df = r2sl.df_rhea2participants2slm.copy()
        # Extract L/R suffix from 'reaction_side' column into new 'side' column.
        df['side'] = df['reaction_side'].str[-1]
        # Determine which rows match a class used in the enumeration.
        df['match'] = df['slm_id'].isin(sl.df_slm_parent2iso['parent_slm_id'])
        # Pivot the dataframe: group by 'MASTER_ID', split by 'side' values, aggregate 'match' values (True/False)
        # -> MASTER_ID | L | R
        # NB: df.pivot_table() is faster than df.apply() -> C‑optimized, no Python loops.
        agg = df.pivot_table(index='MASTER_ID', columns='side', values='match', aggfunc='any', fill_value=False)
        # Get the MASTER_IDs where both sides have at least one 'match'=True.
        l_rhea_used_in_enumeration = agg[(agg['L']) & (agg['R'])].index
    # - Option: At least ONE side has a participant that is a parent class SLMs IDs of isomeric subspecies.
    else:
        l_rhea_used_in_enumeration = r2sl.df_rhea2participants2slm.loc[r2sl.df_rhea2participants2slm['slm_id'].isin(sl.df_slm_parent2iso['parent_slm_id']), 'MASTER_ID'].unique().tolist()
        
    if DEBUG > 1:
        df_rhea_not_used_in_enumeration = df_rhea[~df_rhea['MASTER_ID'].isin(l_rhea_used_in_enumeration)]
        df_rhea_not_used_in_enumeration.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea_not_used_in_enumeration.tsv'), sep="\t", header=True, index=False)
    
    r2sl.df_rhea2participants2slm = r2sl.df_rhea2participants2slm[r2sl.df_rhea2participants2slm['MASTER_ID'].isin(l_rhea_used_in_enumeration)]
    add_statistics(logger, "participants", "%8d Rhea participants (after filtering by parent class SLM ID)\n" % len(r2sl.df_rhea2participants2slm))
    add_statistics(logger, "rhea", "%8d Rhea IDs (after filtering by parent class SLM ID)\n" % len(r2sl.df_rhea2participants2slm['MASTER_ID'].unique()))
    
    #add_statistics(logger, "lipids", "%8d SwissLipids SLMs IDs in Rhea (all levels, before filter_rhea = %s)\n" % (len(SLMs_in_rhea), filter_rhea))
    SLMs_in_rhea_filtered = r2sl.df_rhea2participants2slm['slm_id'].unique()
    #add_statistics(logger, "lipids", "%8d SwissLipids SLMs IDs in Rhea (all levels, after filter_rhea = %s)\n" % (len(SLMs_in_rhea_filtered), filter_rhea))
    
    # Compute and cache the bond changes (we use them later to filter the enumeration results).
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
        # TODO: Do we need this copy?
        df_rhea = df_rhea.copy()
        # TODO: What is this? It seems to be fixing some SMILES?
        df_rhea['rxnsmiles_I'] = df_rhea['rxnsmiles'].apply(lambda x: x.replace('[1*]','C').replace('At','C').replace('[2*]','C').replace('*','C'))
        df_rhea['rxnsmiles_I_mapped'] = df_rhea['rxnsmiles_I'].progress_apply(lambda x: __get_attention_guided_atom_map_error_wrap(x))
        logger.info( "progress_apply: Generating rhea_reactions_bond_changes.tsv - Part 2" )
        df_rhea['bond_changes'] = df_rhea['rxnsmiles_I_mapped'].progress_apply(lambda x: len(r2sl.mapped_reaction_to_report(x)['bond_changes']))
        df_rhea.to_csv(f_cache_rhea_bond_changes, sep='\t', index=False)
        
    logger.info( "Reading %s", f_cache_rhea_bond_changes )
    df_rhea_bond_changes = pd.read_csv(f_cache_rhea_bond_changes, sep='\t')
    logger.info( "LENGTH: %8d df_rhea_bond_changes", len(df_rhea_bond_changes))
    dict_rhea2bond_changes = dict(zip(df_rhea_bond_changes['MASTER_ID'], df_rhea_bond_changes['bond_changes']))

    # Optionally, filter fatty acids.
    df_slm_parent2iso = sl.filter_df_slm_parent2iso(filter_fa, stats, test)

    # This is df_rhea_slm_parent2iso without filters (just for printing statistics).
    df_rhea_slm_parent2iso = sl.df_slm_parent2iso[sl.df_slm_parent2iso['parent_slm_id'].isin(SLMs_in_rhea)]
    add_statistics(logger, "lipids", "-------- Rhea SLMs before filtering\n")
    add_statistics(logger, "lipids", "%8d Rhea pairs of parent class/isomeric subspecies SLM IDs\n" % len(df_rhea_slm_parent2iso))
    add_statistics(logger, "lipids", "%8d Rhea isomeric subspecies SLM IDs\n" % len(set(df_rhea_slm_parent2iso['iso_slm_id'])))
    add_statistics(logger, "lipids", "%8d Rhea parent class SLM IDs\n" % len(set(df_rhea_slm_parent2iso['parent_slm_id'])))

    # This is df_rhea_slm_parent2iso with filter_fa (for next step).
    df_rhea_slm_parent2iso = df_slm_parent2iso[df_slm_parent2iso['parent_slm_id'].isin(SLMs_in_rhea)]
    add_statistics(logger, "lipids", "-------- Rhea SLMs after filter_fa = %s\n" % filter_fa)
    add_statistics(logger, "lipids", "%8d Rhea pairs of parent class/isomeric subspecies SLM IDs\n" % len(df_rhea_slm_parent2iso))
    add_statistics(logger, "lipids", "%8d Rhea isomeric subspecies SLM IDs\n" % len(set(df_rhea_slm_parent2iso['iso_slm_id'])))
    add_statistics(logger, "lipids", "%8d Rhea parent class SLM IDs\n" % len(set(df_rhea_slm_parent2iso['parent_slm_id'])))
    if DEBUG > 1:
        df_rhea_slm_parent2iso.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea_slm_parent2iso.tsv'), sep="\t", header=True, index=False)
    
    # Map Rhea participants parent class to isomeric subspecies SLM IDs.
    r2sl.build_df_rhea2participants2slm_parent2iso(df_rhea_slm_parent2iso, sl.swisslipids)
    add_statistics(logger, "participants", "%8d Rhea participants (after mapping participant parent class SLM ID to isomeric subspecies SLM IDs)\n" % len(r2sl.df_rhea2participants2slm_parent2iso))
    add_statistics(logger, "rhea", "%8d Rhea IDs (after mapping participant parent class SLM ID to isomeric subspecies SLM IDs)\n" % len(r2sl.df_rhea2participants2slm_parent2iso['MASTER_ID'].unique()))

    # This is df_rhea_slm_parent2iso with filter_fa and filer_rhea (just for printing statistics).
    # NB: There are 'nan' values after the LEFT merge, therfore use dropna() to exclude them from the stats.
    df_rhea_slm_parent2iso = (r2sl.df_rhea2participants2slm_parent2iso[['parent_slm_id', 'iso_slm_id']]
                              .dropna(subset=['parent_slm_id', 'iso_slm_id'])
                              .drop_duplicates()
    )
    add_statistics(logger, "lipids", "-------- Rhea SLMs after filter_rhea = %s\n" % filter_rhea)
    add_statistics(logger, "lipids", "%8d Rhea pairs of parent class/isomeric subspecies SLM IDs\n" % len(df_rhea_slm_parent2iso))
    add_statistics(logger, "lipids", "%8d Rhea isomeric subspecies SLM IDs\n" % len(set(df_rhea_slm_parent2iso['iso_slm_id'])))
    add_statistics(logger, "lipids", "%8d Rhea parent class SLM IDs\n" % len(set(df_rhea_slm_parent2iso['parent_slm_id'])))
    if DEBUG > 1:
        df_rhea_slm_parent2iso.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea_slm_parent2iso-after-rhea-filter.tsv'), sep="\t", header=True, index=False)
        
    # Function save_results enumerates reactions and saves the results
    r2sl.save_results(df_rhea2participants, "enumerated_reactions.tsv", dict_rhea2bond_changes, stats)

    # Save statistics to a tsv file.
    save_statistics(output_dir)
