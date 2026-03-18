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
    f_results_statistics = f'{timestamp}_results_statistics.tsv'
    stats = {
        "rhea": [], # unique Rhea IDs
        "reactions": [], # enumerated reactions
        "participants": [],
        "lipids": []
    }
    
    # ---------- Load and Process RheaDB Data ----------

    # Read SwissLipids files.
    sl = SwissLipids(output_dir=output_dir, timestamp=timestamp)
    sl.read_swisslipids_from_file(filter_fa, test)

    # df_slm2chebi maps an SLM ID to a CHEBI ID:
    # Lipid ID | chebi_id
    df_slm2chebi = sl.build_df_slm2chebi()
    stats["lipids"].append( "%8d SwissLipids SLM IDs with a ChEBI ID (df_slm2chebi)\n" % len(df_slm2chebi) )
    if DEBUG > 1:
        logger.debug( "FORMAT: df_slm2chebi\n%s", debug_df_first_row(df_slm2chebi) )
        df_slm2chebi.to_csv(os.path.join(output_dir, 'DEBUG_df_slm2chebi.tsv'), sep="\t", header=True, index=False)

    # Get the unique SLM lipid classes of all isomeric subspecies:
    df_iso_class = sl.get_isomeric_subspecies_class()
    stats["lipids"].append( "%8d SwissLipids class SLM IDs of isomeric subspecies (df_iso_class)\n" % len(df_iso_class) )
    # - Get the ChEBI ID for each SLM ID.
    df_iso_class2chebi = df_iso_class.merge(df_slm2chebi, left_on='Lipid class*', right_on='Lipid ID', how='inner')
    df_iso_class2chebi = df_iso_class2chebi[["Lipid ID", "chebi_id"]]
    stats["lipids"].append( "%8d SwissLipids class ChEBI IDs of isomeric subspecies (df_iso_class2chebi)\n" % len(df_iso_class2chebi) )
    if DEBUG > 1:
        logger.debug( "FORMAT: df_iso_class2chebi\n%s", debug_df_first_row(df_iso_class2chebi) )
        df_iso_class2chebi.to_csv(os.path.join(output_dir, 'DEBUG_df_iso_class2chebi.tsv'), sep="\t", header=True, index=False)
    # - Get the list of ChEBIs.
    l_iso_class_chebi = df_iso_class2chebi['chebi_id'].tolist()
    
    # Read the Rhea FTP files.
    rdb = RheaDB(rhea_version=rhea_version)
    rhea_version = rdb.rhea_db_version
    logger.info( "Loaded files of Rhea release %s", rhea_version )

    # df_rhea has one line per Rhea forward reaction:
    # MASTER_ID | reaction_participant_names | rheaid | rxnsmiles_no_residue_correction | chebi_equation | residue_rxn_flag | rxnsmiles | rxn_smiles_halogen | rxn_smiles_no_A_AH | class_reaction_flag | polymer_reaction_flag | RInChI | Web-RInChIKey
    df_rhea = rdb.df_reactions
    stats["rhea"].append( "%8d Rhea IDs (unfiltered)\n" % len(df_rhea) )
    if DEBUG > 1:
        logger.debug( "FORMAT: df_rhea\n%s", debug_df_first_row(df_rhea) )
        df_rhea.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea.tsv'), sep="\t", header=True, index=False)
    
    # df_rhea2participants has one line per Rhea participant:
    # MASTER_ID | reaction_side | chebiid | smiles | inchi | inchikey | inchikey14L | stoich_coef
    # NB: 'chebiid' also contains POLYMER:x, but for GENERICs it contains the residues' CHEBI IDs.
    df_rhea2participants = rdb.rhea_reaction_long_format_smiles_chebi
    stats["participants"].append( "%8d Rhea participants (unfiltered)\n" % len(df_rhea2participants) )
    if DEBUG > 1:
        logger.debug( "FORMAT: df_rhea2participants\n%s", debug_df_first_row(df_rhea2participants) )
        df_rhea2participants.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea2participants.tsv'), sep="\t", header=True, index=False)
    if rhea_id:
        df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'] == rhea_id]
        
    # Get the subset of Rhea reactions with class compounds (compounds with R-groups).
    # 1. Remove reactions with Web-RInChIKey.
    df_rhea_without_rinchikey = df_rhea[df_rhea['Web-RInChIKey'].isna()]
    df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'].isin(df_rhea_without_rinchikey['MASTER_ID'])]
    stats["participants"].append( "%8d Rhea participants (after filtering by Web-RInChIKey)\n" % len(df_rhea2participants) )
    stats["rhea"].append( "%8d Rhea IDs (after filtering by Web-RInChIKey)\n" % len(df_rhea2participants['MASTER_ID'].unique()) )
    # 2. Remove reactions with Rhea generic compounds
    df_rhea_without_generics  = df_rhea[df_rhea['residue_rxn_flag'] == False]
    df_rhea2participants = df_rhea2participants[df_rhea2participants['MASTER_ID'].isin(df_rhea_without_generics['MASTER_ID'])]
    stats["participants"].append( "%8d Rhea participants (after filtering by generics)\n" % len(df_rhea2participants) )
    stats["rhea"].append( "%8d Rhea IDs (after filtering by generics)\n" % len(df_rhea2participants['MASTER_ID'].unique()) )    
    
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
    stats["participants"].append( "%8d Rhea participants (after filtering by polymers)\n" % len(df_rhea2participants) )
    stats["rhea"].append( "%8d Rhea IDs (after filtering by polymers)\n" % len(df_rhea2participants['MASTER_ID'].unique()) )
    
    # Class RheaToSwissLipids merges Rhea and SwissLipids data.
    r2sl = RheaToSwisslipidsDf(output_dir=output_dir, timestamp=timestamp)
    # df_rhea2participants2slm maps Rhea participants chebi_id to the corresponding SLM ID ('Lipid ID'):
    # MASTER_ID | reaction_side | chebiid | smiles | inchi | inchikey | inchikey14L | stoich_coef | id_prefix | chebi_id | Lipid ID | CHEBI
    r2sl.build_df_rhea2participants2slm(df_rhea2participants, df_slm2chebi)
    stats["participants"].append( "%8d Rhea participants (after mapping participant ChEBI ID to SLM ID)\n" % len(r2sl.df_rhea2participants2slm) )
    stats["rhea"].append( "%8d Rhea IDs (after mapping participant ChEBI ID to SLM ID)\n" % len(r2sl.df_rhea2participants2slm['MASTER_ID'].unique()) )
    if DEBUG > 1:
        logger.debug( "FORMAT: df_rhea2participants2slm\n%s", debug_df_first_row(r2sl.df_rhea2participants2slm) )
        r2sl.df_rhea2participants2slm.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea2participants2slm.tsv'), sep="\t", header=True, index=False)

    # Get the list of unique SwissLipids SLM IDs in Rhea.
    SLMs_in_rhea = r2sl.df_rhea2participants2slm['Lipid ID'].unique()
    stats["lipids"].append( "%8d SwissLipids SLMs IDs in Rhea\n" % len(SLMs_in_rhea) )
    # Check which isomeric subspecies are not covered in Rhea:
    # some are candidates for new reactions?
    df_iso_class_in_rhea     = df_iso_class2chebi[df_iso_class2chebi['Lipid ID'].isin(SLMs_in_rhea)][["Lipid ID", "chebi_id"]]
    df_iso_class_not_in_rhea = df_iso_class2chebi[~df_iso_class2chebi['Lipid ID'].isin(SLMs_in_rhea)][["Lipid ID", "chebi_id"]]
    stats["lipids"].append( "%8d SwissLipids class SLMs IDs of isomeric subspecies in Rhea\n"     % len(df_iso_class_in_rhea) )
    stats["lipids"].append( "%8d SwissLipids class SLMs IDs of isomeric subspecies not in Rhea\n" % len(df_iso_class_not_in_rhea) )
    if DEBUG > 1:
        df_iso_class_in_rhea.to_csv(os.path.join(output_dir, 'DEBUG_df_iso_class_in_rhea.tsv'), sep="\t", header=True, index=False) 
        df_iso_class_not_in_rhea.to_csv(os.path.join(output_dir, 'DEBUG_df_iso_class_not_in_rhea.tsv'), sep="\t", header=True, index=False)
    
    # Keep only reactions with at least one participant that was used as a class in the enumeration.
    # NB: We can only do this here, because it needs a clean column 'chebi_id'.
    if filter_rhea == True:
        l_rhea_used_in_enumeration = r2sl.df_rhea2participants2slm.loc[r2sl.df_rhea2participants2slm['chebi_id'].isin(l_iso_class_chebi), 'MASTER_ID'].unique().tolist()
        if DEBUG > 1:
            df_rhea_not_used_in_enumeration = df_rhea[~df_rhea['MASTER_ID'].isin(l_rhea_used_in_enumeration)]
            df_rhea_not_used_in_enumeration.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea_not_used_in_enumeration.tsv'), sep="\t", header=True, index=False)
        
        r2sl.df_rhea2participants2slm = r2sl.df_rhea2participants2slm[r2sl.df_rhea2participants2slm['MASTER_ID'].isin(l_rhea_used_in_enumeration)]
        stats["participants"].append( "%8d Rhea participants (after filtering by SLM used in enumeration)\n" % len(r2sl.df_rhea2participants2slm) )
        stats["rhea"].append( "%8d Rhea IDs (after filtering by SLM used in enumeration)\n" % len(r2sl.df_rhea2participants2slm['MASTER_ID'].unique()) )
        SLMs_in_rhea = r2sl.df_rhea2participants2slm['Lipid ID'].unique()
        stats["lipids"].append( "%8d SwissLipids SLMs IDs in Rhea (after filtering by SLM used in enumeration)\n" % len(SLMs_in_rhea) )    
        
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
    logger.info( "LENGTH: %8d df_rhea_bond_changes", len(df_rhea_bond_changes) )
    dict_rhea2bond_changes = dict(zip(df_rhea_bond_changes['MASTER_ID'], df_rhea_bond_changes['bond_changes']))

    # df_rhea_slm_class2iso maps a Rhea class SLM ID to its isomeric subspecies SLM IDs:
    # class_slm_id | iso_slm_id | Lipid ID | Level | Name | Lipid class* | Components* | SMILES (pH7.3) | CHEBI | sn1' | sn2' | ...
    # NB: iso_slm_id = Lipid ID        
    # Optionally, filter fatty acids.
    df_rhea_slm_class2iso = sl.build_and_filter_df_slm_class2iso(SLMs_in_rhea, filter_fa, stats, test)
    stats["lipids"].append( "-------- df_rhea_slm_class2iso statistics:\n" )
    stats["lipids"].append( "%8d pairs of Rhea class/isomeric subspecies SLM IDs\n" % len(df_rhea_slm_class2iso) )
    stats["lipids"].append( "%8d Rhea class SLM IDs\n" % len(set(df_rhea_slm_class2iso["class_slm_id"])) )
    stats["lipids"].append( "%8d Rhea isomeric subspecies SLM IDs\n" % len(set(df_rhea_slm_class2iso["iso_slm_id"])) )
    if DEBUG > 1:
        logger.debug( "FORMAT: df_rhea_slm_class2iso\n%s", debug_df_first_row(df_rhea_slm_class2iso) )
        df_rhea_slm_class2iso.to_csv(os.path.join(output_dir, 'DEBUG_df_rhea_slm_class2iso.tsv'), sep="\t", header=True, index=False)
    
    # Map Rhea participants class to isomeric subspecies SLM IDs.
    r2sl.build_df_rhea2participants2slm_class2iso(df_rhea_slm_class2iso, sl.swisslipids)
    stats["participants"].append( "%8d Rhea participants (after mapping participant class SLM ID to isomeric subspecies SLM IDs)\n" % len(r2sl.df_rhea2participants2slm_class2iso) )
    stats["rhea"].append( "%8d Rhea IDs (after mapping participant class SLM ID to isomeric subspecies SLM IDs)\n" % len(r2sl.df_rhea2participants2slm_class2iso['MASTER_ID'].unique()) )
    # NB: There are 'nan' values after the LEFT merge, therfore use dropna() to exclude them from the stats.
    stats["lipids"].append( "%8d Rhea class SLM IDs (after mapping participant class SLM ID to isomeric subspecies SLM IDs)\n" % len(set(r2sl.df_rhea2participants2slm_class2iso['class_slm_id'].dropna())) )
    stats["lipids"].append( "%8d Rhea isomeric subspecies SLM IDs (after mapping participant class SLM ID to isomeric subspecies SLM IDs)\n" % len(set(r2sl.df_rhea2participants2slm_class2iso['iso_slm_id'].dropna())) )
        
    # Function save_results generated the combinations and saves the results
    r2sl.save_results(df_rhea2participants, f'{timestamp}_enumerated_reactions.tsv', dict_rhea2bond_changes, stats)

    # Save statistics to a tsv file.
    with open(os.path.join(output_dir, f_results_statistics), 'w') as f:
        
        f.write("\n\nRhea ID statistics:\n\n")
        for line in stats["rhea"]:
            f.write(line)

        f.write("\n\nRhea participants statistics:\n\n")
        for line in stats["participants"]:
            f.write(line)

        f.write("\n\nSLM lipids statistics:\n\n")
        for line in stats["lipids"]:
            f.write(line)

        f.write("\n\nEnumerated reactions statistics:\n\n")
        for line in stats["reactions"]:
            f.write(line)

    logger.info( "Statistics saved to %s\n", f_results_statistics )
