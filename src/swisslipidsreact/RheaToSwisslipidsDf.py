import logging
import os
import re
import signal
import pandas as pd
import numpy as np
from itertools import product
from collections import Counter
from typing import NamedTuple, List, Dict, Iterable, Tuple, Optional
from tqdm import tqdm

from rxnmapper import RXNMapper
from rdkit import Chem

from pyrheadb.Reaction import Reaction
from pyrheadb.RInChI import RInChI

from .utils import DEBUG, debug_df_first_row
from .FA_lists import positions

logger = logging.getLogger(__name__)

# tqdm.pandas() wraps pandas.apply() method with a tqdm progress bar.
# Use for long-running apply() if you want to monitor progress.
tqdm.pandas()

class TimeoutException(Exception):
    pass

def handler(signum, frame):
    raise TimeoutException()

# Set the signal handler
signal.signal(signal.SIGALRM, handler)

# Lightweight helper class to avoid pandas overhead.
class Participant(NamedTuple):
    chebi_id: str                # row['chebi_id'] = grouping key
    chebi_display: str           # row['CHEBI'] if present else row['chebi_id']
    side: str                    # 'L' or 'R'
    stoich: int                  # int(row['stoich_coef']) or 1
    smiles: str                  # row['SMILES (pH7.3)'] if present else row['smiles']
    iso_lipid_id: Optional[str]  # row['isomeric_subspecies_descendant_lipid_id'] or None
    components: Tuple[str, ...]  # tuple of components; duplicates already included by stoich
    comp_counter: Counter        # Counter of components, e.g. Counter({'C1': 2, 'C2': 1}), for equality check.
    comp_text: str               # str(row['Components*']) for components_equation

class RheaToSwisslipidsDf():
    
    def __init__(self, output_dir='', timestamp='now'):
        self.rxn_mapper = RXNMapper()
        self.output_dir = output_dir
        self.timestamp  = timestamp

    def get_df_swiss_lipids_chebi_rhea(self, df_swiss_lipids_chebi, df_rhea_chebi):
        self.df_swiss_lipids_chebi_rhea = df_rhea_chebi.merge(df_swiss_lipids_chebi, on='chebi_id', how='inner')

    # --- Merge to obtain Rhea reactions of isomeric subspecies ---
    def get_df_rhea_descendant(self, rhea_lipid_to_descendant_df, df_swisslipids):

        # FORMAT: rhea_lipid_to_descendant_df = same as in main.py
        logger.info(  "LENGTH: rhea_lipid_to_descendant_df: %d", len(rhea_lipid_to_descendant_df) )
        logger.info(  "LENGTH: df_swisslipids: %d", len(df_swisslipids) )
        logger.debug( "FORMAT: df_swisslipids\n%s", debug_df_first_row(df_swisslipids) )

        logger.debug( "LEFT merging df_swiss_lipids_chebi_rhea WITH rhea_lipid_to_descendant_df ON 'Lipid ID' = 'rhea_lipid_id'" )
        df_rhea_descendant = self.df_swiss_lipids_chebi_rhea.merge(
            rhea_lipid_to_descendant_df, left_on='Lipid ID', right_on='rhea_lipid_id', how='left'
        )
        logger.info(  "LENGTH: df_rhea_descendant: %d", len(df_rhea_descendant) )
        logger.debug( "FORMAT: df_rhea_descendant\n%s", debug_df_first_row(df_rhea_descendant) )

        logger.debug( "Selecting a subset of columns of df_rhea_descendant" )
        df_rhea_descendant = df_rhea_descendant[[
            'MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef',
            'rhea_lipid_id', 'isomeric_subspecies_descendant_lipid_id', 'smiles'
        ]]
        logger.debug( "FORMAT: df_rhea_descendant\n%s", debug_df_first_row(df_rhea_descendant) )

        logger.debug( "LEFT merging df_rhea_descendant WITH df_swisslipids ON 'Lipid ID'" )
        df_rhea_descendant = df_rhea_descendant.merge(
            df_swisslipids, left_on='isomeric_subspecies_descendant_lipid_id',
            right_on='Lipid ID', how='left'
        )
        logger.info(  "LENGTH: df_rhea_descendant: %d", len(df_rhea_descendant) )
        logger.debug( "FORMAT: df_rhea_descendant\n%s", debug_df_first_row(df_rhea_descendant) )
        
        logger.debug( "Dropping 3 columns of df_rhea_descendant: 'Lipid ID', 'Level', 'Lipid class*'" )
        df_rhea_descendant.drop(columns=['Lipid ID', 'Level', 'Lipid class*'], inplace=True)
        logger.info(  "LENGTH: df_rhea_descendant: %d", len(df_rhea_descendant) )
        logger.debug( "FORMAT: df_rhea_descendant\n%s", debug_df_first_row(df_rhea_descendant) )

        self.df_rhea_descendant = df_rhea_descendant
        return len(df_rhea_descendant), len(set(df_rhea_descendant['isomeric_subspecies_descendant_lipid_id']))


    # --- Reaction generation functions ---
    def __build_equations(self, left_combination: Iterable[Participant], right_groups: List[List[Participant]]):
        """
        Constructs reaction SMILES and equations for a given left-side combination
        against all right-side combinations. Pure-Python, no per-iteration DataFrames.
        Returns (list_of_4tuples, 'yes_component_match' | 'no_component_match').
        """
        # Helper function to build equations from left and right list of Participants.
        def __build_equations_from_participants(left_participants: List[Participant], right_participants: List[Participant]) -> Tuple[str, str, str, str]:
            smiles_left: List[str] = []
            smiles_right: List[str] = []
            chebi_left: List[str] = []
            chebi_right: List[str] = []
            sl_left: List[str] = []
            sl_right: List[str] = []
            comp_left: List[str] = []
            comp_right: List[str] = []

            # Process both sides in one pass to minimize Python overhead
            for row in left_participants:
                if row.stoich > 0 and row.smiles:
                    smiles_left.extend([row.smiles] * row.stoich)
                chebi_left.append(f"{row.stoich} {row.chebi_display}")
                sl_left.append(f"{row.stoich} {row.iso_lipid_id if row.iso_lipid_id is not None else 'NA'}")
                comp_left.append(row.comp_text)

            for row in right_participants:
                if row.stoich > 0 and row.smiles:
                    smiles_right.extend([row.smiles] * row.stoich)
                chebi_right.append(f"{row.stoich} {row.chebi_display}")
                sl_right.append(f"{row.stoich} {row.iso_lipid_id if row.iso_lipid_id is not None else 'NA'}")
                comp_right.append(row.comp_text)

            reaction_smiles      = '.'.join(smiles_left) + '>>'  + '.'.join(smiles_right)
            chebi_equation      = ' + '.join(chebi_left) + ' = ' + ' + '.join(chebi_right)
            sl_equation         = ' + '.join(sl_left)    + ' = ' + ' + '.join(sl_right)
            components_equation = ' + '.join(comp_left)  + ' = ' + ' + '.join(comp_right)
            return reaction_smiles, chebi_equation, sl_equation, components_equation

        # --- Main section ---
        res_final:  List[Tuple[str, str, str, str]] = []
        res_backup: List[Tuple[str, str, str, str]] = []

        # Precompute the frequency of all components of the left participants (as a multiset).
        left_counter = Counter()
        left_components: set = set()
        left_participants: List[Participant] = list(left_combination)
        for participant in left_participants:
            left_counter += participant.comp_counter
            if participant.components:
                left_components.update(participant.components)

        # Filter the candidates of right participants by the left participants' components (if any).
        right_groups_filtered: List[List[Participant]] = []
        for right_participants in right_groups:
            if left_components:
                filtered_group = [participant
                                  for participant in right_participants
                                  if (not participant.components)
                                  or all((comp in left_components) for comp in participant.components)]
            else:
                filtered_group = right_participants

            # If the filtered right side is empty, return empty list + no match
            if not filtered_group:
                return [], 'no_component_match'
            
            right_groups_filtered.append(filtered_group)
        
        # Enumerate the right-side combinations and eliminate those with unequal components counts.
        for right_combination in product(*right_groups_filtered):
            right_counter = Counter()
            for participant in right_combination:
                right_counter += participant.comp_counter

            # Get the total number of components on each side.    
            left_size  = sum(left_counter.values())
            right_size = sum(right_counter.values())

            # Good combinations are those with components on the left side
            # and where the number of components are equal on both sides.
            if left_size > 0 and left_counter == right_counter:
                res_final.append(__build_equations_from_participants(left_participants, list(right_combination)))
            else: 
                res_backup.append(__build_equations_from_participants(left_participants, list(right_combination)))

        if res_final:
            return res_final, 'yes_component_match'
        return res_backup, 'no_component_match'

    def __enumerate_reactions(self, df):
        """Generates all valid reaction SMILES and equation strings for each MASTER_ID."""

        # Helper function to get stoichiometric list of components.
        def __extract_stoichiometric_components(row):
            """
            Builds an array of all components, with each component repeated as many times
            as the stoichiometric coefficient of the reaction participant.
            """
            stoich = int(row['stoich_coef']) if pd.notnull(row['stoich_coef']) else 1 # Default to 1 if missing.
            components = []
            # Loop over the columns where the components are stored.
            for col in positions:
                if col in row and pd.notnull(row[col]) and isinstance(row[col], str) and row[col].strip():
                    components.extend([row[col]] * stoich)
            return components
        
        # Helper function to generate a Participant object.
        def __get_participant(rec: Dict) -> Participant:
            chebi_id = rec['chebi_id']
            chebi_display = rec.get('CHEBI')
            if chebi_display is None or (isinstance(chebi_display, float) and pd.isna(chebi_display)):
                chebi_display = rec['chebi_id']
            side = rec['side'] # 'L' or 'R'
            stoich = int(rec.get('stoich_coef', 1) or 1)
            smiles_pH = rec.get('SMILES (pH7.3)')
            smiles = smiles_pH if (smiles_pH is not None and pd.notna(smiles_pH)) else rec.get('smiles', '')
            iso_lipid_id = rec.get('isomeric_subspecies_descendant_lipid_id')
            components = tuple(rec.get('components', []) or [])
            comp_counter = rec.get('comp_counter') or Counter(components)
            comp_text = str(rec.get('Components*'))
            return Participant(chebi_id, chebi_display, side, stoich, smiles, iso_lipid_id, components, comp_counter, comp_text)

        # Helper function to group participants by their class.
        def __group_participants_by_class(df_participants: pd.DataFrame) -> List[List[Participant]]:
            """
            Groups the participants (isomeric subspecies lipid ID) by their class (CHEBI ID).
            Returns list of participants lists, grouped by CHEBI ID, keeping stable order.
            """
            groups: Dict[str, List[Participant]] = {}
            # Convert DataFrame to Dict of lipid class(CHEBI) -> list of isomeric subspecies(participants).
            for rec in df_participants.to_dict('records'):
                participant = __get_participant(rec)
                groups.setdefault(participant.chebi_id, []).append(participant)
            # Keep insertion order of first occurrence of CHEBI without sorting.
            return list(groups.values())

        # --- Main section ---
        logger.info("progress_apply: Extracting stoichiometric coefficients per component")
        # Get an array of all components, where each component appears as many times as its stoichiometric coefficient.
        df['components'] = df.progress_apply(__extract_stoichiometric_components, axis=1)
        # Precompute how many times each component occurs
        # (we will use this for pruning impossible combinations).
        df['comp_counter'] = df['components'].apply(Counter)
        # Precompute side.
        df['side'] = df['reaction_side'].apply(lambda x: x.split('_')[1])
        df_left  = df[df['side'] == 'L']
        df_right = df[df['side'] == 'R']

        total_gen_attempted = len(df_left['reaction_side'].unique())
        results = []
        
        logger.info("progress_apply: Enumerating reactions")
        for reaction_side, df_left_group in tqdm(df_left.groupby('reaction_side')):
            master_id = int(reaction_side.split('_')[0])
            df_right_group = df_right[df_right['MASTER_ID'] == master_id]

            signal.alarm(90000)
            try:
                # Convert pandas to lightweight Python data structure to make the enumeration fast.  
                left_groups:  List[List[Participant]] = __group_participants_by_class(df_left_group)
                right_groups: List[List[Participant]] = __group_participants_by_class(df_right_group)

                # If any side has no groups, there is nothing to enumerate.
                if not left_groups or not right_groups:
                    continue

                # Loop over the cartesian product of the left combinations:
                for left_combination in product(*left_groups):
                    # For each left combination, build the equations with the right combinations with pruning.
                    res_total, component_match = self.__build_equations(left_combination, right_groups)
                    for reaction_smiles, chebi_equation, sl_equation, components_equation in res_total:
                        if reaction_smiles:
                            results.append({
                                'MASTER_ID': master_id,
                                'chebi_equation': chebi_equation,
                                'swisslipids_equation': sl_equation,
                                'components_equation': components_equation,
                                'reaction_smiles': reaction_smiles,
                                'component_match': component_match
                            })
            except TimeoutException:
                logger.warning("MASTER ID %s took too long and was skipped: %d",
                               reaction_side_master_id.split('_')[0], len(df))
                continue
            finally:
                signal.alarm(0) # Cancel alarm.

        return pd.DataFrame(results), total_gen_attempted


    def mapped_reaction_to_report(self, mapped_rxn):

        def __extract_atom_environment(mol):
            """Return a dict of atom map number -> (symbol, list of (neighbor map num, bond type))"""
            env = {}
            for atom in mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num == 0:
                    continue
                neighbors = []
                for bond in atom.GetBonds():
                    nbr = bond.GetOtherAtom(atom)
                    nbr_map = nbr.GetAtomMapNum()
                    if nbr_map > 0:
                        neighbors.append((nbr_map, str(bond.GetBondType())))
                env[map_num] = {
                    'symbol': atom.GetSymbol(),
                    'neighbors': sorted(neighbors)
                }
            return env      
        
        reactant_smiles, product_smiles = mapped_rxn.split(">>")
        reactants = [Chem.MolFromSmiles(s) for s in reactant_smiles.split(".")]
        products  = [Chem.MolFromSmiles(s) for s in product_smiles.split(".")]

        react_env = {}
        for mol in reactants:
            react_env.update(__extract_atom_environment(mol))

        prod_env = {}
        for mol in products:
            prod_env.update(__extract_atom_environment(mol))

        all_map_nums = set(react_env.keys()).union(prod_env.keys())

        report = {
            'missing_in_products': [],
            'missing_in_reactants': [],
            'element_mismatches': [],
            'bond_changes': [],
            'num_preserved_atoms': 0,
            'num_total_mapped_atoms': len(all_map_nums)
        }

        for num in sorted(all_map_nums):
            r = react_env.get(num)
            p = prod_env.get(num)

            if r and not p:
                report['missing_in_products'].append((num, r['symbol']))
            elif p and not r:
                report['missing_in_reactants'].append((num, p['symbol']))
            elif r['symbol'] != p['symbol']:
                report['element_mismatches'].append((num, r['symbol'], p['symbol']))
            else:
                report['num_preserved_atoms'] += 1
                if r['neighbors'] != p['neighbors']:
                    report['bond_changes'].append({
                        'map_num': num,
                        'symbol': r['symbol'],
                        'reactant_bonds': r['neighbors'],
                        'product_bonds': p['neighbors']
                    })

        return report

    # ---------- Generate, Filter, and Save Final Results ----------
    def save_results(self, rheadf, df, rhea_reactions, filename, bond_changes_lookup):

        # Helper function
        def __get_position2component_map(components_equation_side: str) -> dict[str, str]:
            """
            Parses one side of a components_equation string to make a map of position->component.
            Uses the same pattern and parsing logic as SwissLipids.py.
            """
            # pattern like "SLM:000000123 (sn1 or sn2)"
            pattern = r'(SLM:\d+)\s+\(([^)]+)\)'
            matches = re.findall(pattern, components_equation_side)

            extracted_rows = []
            for slm_id, pos_string in matches:
                pos_split = re.split(r'\s*(?:or|and|,)\s*', pos_string)
                for slm_pos in pos_split:
                    slm_pos = slm_pos.strip()
                    if slm_pos in positions:
                        extracted_rows.append((slm_pos, slm_id))

            # Keep only one SLM (the first) per position.
            pos_map = {}
            for slm_pos, slm_id in extracted_rows:
                if slm_pos not in pos_map:
                    pos_map[slm_pos] = slm_id
            return pos_map

        # Helper function
        def __get_component_changes(components_equation: str) -> dict[str, str]:
            """
            Parses a components_equation string and compares for each position
            the component on the left- and right-hand side and stores whether
            the component has changed.

            Example: SLM:000000510 (sn2) / SLM:000001101 (sn1) = nan + SLM:000000510 (sn1 or sn2)
                     sn1 : SLM:000001101-SLM:000000510
                     sn1_change : True
                     sn2 : SLM:000000510-SLM:000000510
                     sn2_change : False
            """
            left_string, right_string = components_equation.split('=', 1)
            left_map  = __get_position2component_map(left_string)
            right_map = __get_position2component_map(right_string)

            result = {}
            for p in positions:
                left  = left_map.get(p)
                right = right_map.get(p)

                if left is None and right is None:
                    result[p] = ''
                    result[f'{p}_change'] = False
                else:
                    result[p] = f'{left or "None"}-{right or "None"}'
                    # If only one side has a component, we consider this OK an set {p}_change = False.
                    result[f'{p}_change'] = (left is not None and right is not None and left != right)
            return result

        # Helper function
        def __summarize_discrepancies(reaction_smiles):
            try:
                result = self.rxn_mapper.get_attention_guided_atom_maps([reaction_smiles])[0]
            except:
                report = {
                    'missing_in_products': [],
                    'missing_in_reactants': [],
                    'element_mismatches': [],
                    'bond_changes': [],
                    'num_preserved_atoms': 0,
                    'num_total_mapped_atoms': 0
                }
                return report
            mapped_rxn = result['mapped_rxn']
            report = self.mapped_reaction_to_report(mapped_rxn)
            return report

        # Helper function
        def __check_reactions_for_which_component_matching_was_impossible_with_atom_mapper(row, bond_changes_lookup):
            if row['component_match'] == 'no_component_match':
                discrepancies = __summarize_discrepancies(row['reaction_smiles'])
                return len(discrepancies['bond_changes']) <= bond_changes_lookup[row['MASTER_ID']]
            return True
        

        # --- Main section ---
        # Arguments: df = r2sl.df_rhea_descendant
        logger.info(  "LENGTH: rheadf: %d", len(rheadf) )
        logger.info(  "LENGTH: df: %d", len(df) )
        logger.debug( "FORMAT: df\n%s", debug_df_first_row(df) )
        if DEBUG > 1:
            df.to_csv(os.path.join(self.output_dir, 'DEBUG_df.tsv'), sep="\t", header=True, index=False)
        
        # Get subset of columns of rhea2participant dataframe.
        df_rhea_cut = rheadf[['MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'smiles']]
        # Get subset of rows where MASTER_ID is in df (= r2sl.df_rhea_descendant).
        df_rhea_cut = df_rhea_cut[df_rhea_cut['MASTER_ID'].isin(df['MASTER_ID'])]
        logger.info(  "LENGTH: df_rhea_cut: %d", len(df_rhea_cut) )
        logger.debug( "FORMAT: df_rhea_cut\n%s", debug_df_first_row(df_rhea_cut) )

        # TODO: Why RIGHT merge? Why first filter df_rhea_cut by MASTER_ID?
        logger.debug( "RIGHT merging df WITH df_rhea_cut ON MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'smiles'" )
        df = df.merge(
            df_rhea_cut, on=['MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'smiles'], how='right'
        )
        logger.info(  "LENGTH: df: %d", len(df) )
        logger.debug( "FORMAT: df\n%s", debug_df_first_row(df) )
        if DEBUG > 0:
            df.to_csv(os.path.join(self.output_dir, 'DEBUG_df-after_merge.tsv'), sep="\t", header=True, index=False)

        result_df, total_gen_attempted = self.__enumerate_reactions(df)
        logger.info( "Statistics: Reactions enumerated: %d\n", len(result_df) )
        logger.debug( "FORMAT: result_df\n%s", debug_df_first_row(result_df) )
        if DEBUG > 1:
            result_df.to_csv(os.path.join(self.output_dir, 'DEBUG_result_df-reactions_enumerated.tsv'), sep="\t", header=True, index=False)
        rhea_reactions['5_after_enumerating_reactions'] = rhea_reactions['MASTER_ID'].isin(result_df['MASTER_ID'])

        # Drop rows without reaction SMILES.
        result_df.dropna(subset = ['reaction_smiles'], inplace=True)
        logger.info( "Statistics: Reactions after dropping those without reaction SMILES: %d\n", len(result_df) )
        rhea_reactions['6_after_dropping_na_reaction_smiles'] = rhea_reactions['MASTER_ID'].isin(result_df['MASTER_ID'])
        MASTER_ID_for_specific_fatty_acids_before_filtering_out_the_unbalanced = len(set(result_df['MASTER_ID']))
        num_reactions_to_check_for_balance = len(result_df)
        
        logger.info( "progress_apply: Checking components balanced" )
        # Example input rows:
        # MASTER_ID components_equation
        # 10100 nan + SLM:000000510 (acyl) = nan + nan + SLM:000000510 (free fatty acid)
        # 10272 SLM:000000510 (sn1 or sn2) + nan = nan + SLM:000000510 (sn1 or sn2) + nan
        # 10352 nan + nan = nan + nan
        parsed = result_df['components_equation'].progress_apply(__get_component_changes).apply(pd.Series)
        result_df = pd.concat([result_df, parsed], axis=1)
        if DEBUG > 0:
            result_df.to_csv(os.path.join(self.output_dir, 'DEBUG_result_df-before-applying-components-balanced.tsv'), sep="\t", header=True, index=False)

        # Store exceptions:
        # 63596, 78043, 77759, 78435 have triglyceride twice in reactants, therefore sn1 is not unique, 
        # there is sn1 of the first tryglyceride, and sn2 of the second, which makes the code confused
        df_change_of_sn_by_design = result_df[result_df['MASTER_ID'].isin([63596, 78043, 77759, 78435])]
        df_change_of_sn_by_design['component_match']='no_component_match'
        # Get list of boolean '_change' columns and remove rows where any of them is True.
        change_cols = [col for col in result_df.columns if col.endswith('_change')]
        result_df = result_df[~result_df[change_cols].any(axis=1)]
        # Now add back the exceptions.
        result_df = pd.concat([result_df, df_change_of_sn_by_design])
        if DEBUG > 0:
            result_df.to_csv(os.path.join(self.output_dir, 'DEBUG_result_df-after-checking-components-balanced.tsv'), sep="\t", header=True, index=False)
        logger.info( "Statistics: Reactions after filtering by component match: %d\n", len(result_df) )
        rhea_reactions['7_after_dropping_incorrect_component_matches'] = rhea_reactions['MASTER_ID'].isin(result_df['MASTER_ID'])

        # Filter for balanced reactions
        rxn = Reaction()
        logger.info( "progress_apply: Checking reaction balanced" )
        result_df['balanced'] = result_df['reaction_smiles'].progress_apply(rxn.check_reaction_balance)
        result_df = result_df[result_df['balanced']==True]
        result_df.drop(columns=['balanced'], inplace=True)
        MASTER_ID_for_specific_fatty_acids_after_filtering_out_the_unbalanced = len(set(result_df['MASTER_ID']))
        rhea_reactions['8_after_dropping_unbalanced_reactions'] = rhea_reactions['MASTER_ID'].isin(result_df['MASTER_ID'])
        num_reactions_after_filtering_out_the_unbalanced = len(result_df)
        logger.info( "Statistics: Reactions after filtering by balance: %d\n", len(result_df) )
        # Filter by bond changes
        logger.info( "progress_apply: Checking number of bond changes" )
        result_df['bond_breakage_max_4'] = result_df.progress_apply(__check_reactions_for_which_component_matching_was_impossible_with_atom_mapper, axis=1, args=[bond_changes_lookup,])
        result_df = result_df[result_df['bond_breakage_max_4']==True]
        logger.info( "Statistics: Reactions after filtering by bond changes: %d\n", len(result_df) )
        rhea_reactions.to_csv(os.path.join(self.output_dir, f'{self.timestamp}_rhea_reactions_overview.tsv'), sep='\t', index=False)

        # Generate RInChI and Web-RInChIKey
        logger.info( "progress_apply: Generating RInChI" )
        rinchi = RInChI()
        if len(result_df)>0:
            result_df[['RInChI', 'Web-RInChIKey']] = result_df.progress_apply(
                lambda x: rinchi.error_handle_wrap_rinchi(x['reaction_smiles']),
                axis=1, result_type='expand'
            )
        
            # Remove reactions without Web-RInChIKey (without structures) and duplicates
            result_df.dropna(subset=['Web-RInChIKey'], inplace=True)
            result_df.drop_duplicates(subset=['MASTER_ID', 'Web-RInChIKey', 'reaction_smiles'], inplace=True)
            logger.info( "Statistics: Number of reactions after removing empty and duplicated by Web-RInChIKey: %d", len(result_df) )
            
            # Save final result
            result_df.to_csv(os.path.join(self.output_dir, filename), sep='\t', index=False)

        return total_gen_attempted, \
              MASTER_ID_for_specific_fatty_acids_before_filtering_out_the_unbalanced, \
               num_reactions_to_check_for_balance, \
               MASTER_ID_for_specific_fatty_acids_after_filtering_out_the_unbalanced, \
               num_reactions_after_filtering_out_the_unbalanced
