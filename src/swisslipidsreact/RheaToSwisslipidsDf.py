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
    iso_slm_id: Optional[str]    # row['iso_slm_id'] or row['chebi_display']
    iso_slm_name: Optional[str]  # row['iso_slm_name'] or row['chebi_display']
    components: Tuple[str, ...]  # tuple of components; duplicates already included by stoich
    comp_counter: Counter        # Counter of components, e.g. Counter({'C1': 2, 'C2': 1}), to prune combinations.
    comp_text: str               # str(row['Components*']) for components_equation

class RheaToSwisslipidsDf():
    
    def __init__(self, output_dir='', timestamp='now'):
        self.rxn_mapper = RXNMapper()
        self.output_dir = output_dir
        self.timestamp  = timestamp

    def build_df_rhea2participants2slm(self, df_rhea2participants, df_slm2chebi):
        self.df_rhea2participants2slm = df_rhea2participants.merge(df_slm2chebi, on='chebi_id', how='inner')

    # --- Merge to obtain Rhea reactions of isomeric subspecies ---
    def build_df_rhea2participants2slm_class2iso(self, df_rhea_slm_class2iso, df_swisslipids):

        # FORMAT: df_rhea_slm_class2iso = same as in main.py
        logger.info( "LENGTH: %8d df_rhea_slm_class2iso", len(df_rhea_slm_class2iso) )
        logger.info( "LENGTH: %8d df_swisslipids", len(df_swisslipids) )
        logger.debug( "FORMAT: df_swisslipids\n%s", debug_df_first_row(df_swisslipids) )

        logger.debug( "LEFT merging self.df_rhea2participants2slm WITH df_rhea_slm_class2iso ON 'Lipid ID' = 'class_slm_id'" )
        df_rhea2participants2slm_class2iso = self.df_rhea2participants2slm.merge(
            df_rhea_slm_class2iso, left_on='Lipid ID', right_on='class_slm_id', how='left'
        )
        logger.info( "LENGTH: %8d df_rhea2participants2slm_class2iso", len(df_rhea2participants2slm_class2iso) )
        logger.debug( "FORMAT: df_rhea2participants2slm_class2iso\n%s", debug_df_first_row(df_rhea2participants2slm_class2iso) )

        logger.debug( "Selecting a subset of columns of df_rhea2participants2slm_class2iso and renaming some" )
        df_rhea2participants2slm_class2iso = df_rhea2participants2slm_class2iso[[
            'MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'class_slm_id', 'iso_slm_id', 'smiles', 'Name'
        ]].rename(columns={'Name': 'iso_slm_name'})
        logger.debug( "FORMAT: df_rhea2participants2slm_class2iso\n%s", debug_df_first_row(df_rhea2participants2slm_class2iso) )

        logger.debug( "LEFT merging df_rhea2participants2slm_class2iso WITH df_swisslipids ON 'iso_slm_id' = 'Lipid ID'" )
        df_rhea2participants2slm_class2iso = df_rhea2participants2slm_class2iso.merge(
            df_swisslipids, left_on='iso_slm_id', right_on='Lipid ID', how='left'
        )
        logger.info( "LENGTH: %8d df_rhea2participants2slm_class2iso", len(df_rhea2participants2slm_class2iso) )
        logger.debug( "FORMAT: df_rhea2participants2slm_class2iso\n%s", debug_df_first_row(df_rhea2participants2slm_class2iso) )
        
        logger.debug( "Dropping 3 columns of df_rhea2participants2slm_class2iso: 'Lipid ID', 'Level', 'Lipid class*'" )
        df_rhea2participants2slm_class2iso.drop(columns=['Lipid ID', 'Level', 'Lipid class*'], inplace=True)
        logger.debug( "FORMAT: df_rhea2participants2slm_class2iso\n%s", debug_df_first_row(df_rhea2participants2slm_class2iso) )

        self.df_rhea2participants2slm_class2iso = df_rhea2participants2slm_class2iso
 

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
            slm_id_left: List[str] = []
            slm_id_right: List[str] = []
            slm_name_left: List[str] = []
            slm_name_right: List[str] = []
            comp_left: List[str] = []
            comp_right: List[str] = []

            # Process both sides in one pass to minimize Python overhead
            for row in left_participants:
                if row.stoich > 0 and row.smiles:
                    smiles_left.extend([row.smiles] * row.stoich)
                chebi_left.append(f"{row.stoich} {row.chebi_display}")
                slm_id_left.append(f"{row.stoich} {row.iso_slm_id if pd.notna(row.iso_slm_id) else row.chebi_display}")
                slm_name_left.append(f"{row.stoich} {row.iso_slm_name if pd.notna(row.iso_slm_name) else row.chebi_display}")
                comp_left.append(row.comp_text)

            for row in right_participants:
                if row.stoich > 0 and row.smiles:
                    smiles_right.extend([row.smiles] * row.stoich)
                chebi_right.append(f"{row.stoich} {row.chebi_display}")
                slm_id_right.append(f"{row.stoich} {row.iso_slm_id if pd.notna(row.iso_slm_id) else row.chebi_display}")
                slm_name_right.append(f"{row.stoich} {row.iso_slm_name if pd.notna(row.iso_slm_name) else row.chebi_display}")
                comp_right.append(row.comp_text)

            reaction_smiles     = '.'.join(smiles_left)     + '>>'  + '.'.join(smiles_right)
            chebi_equation      = ' + '.join(chebi_left)    + ' = ' + ' + '.join(chebi_right)
            slm_id_equation     = ' + '.join(slm_id_left)   + ' = ' + ' + '.join(slm_id_right)
            slm_name_equation   = ' + '.join(slm_name_left) + ' = ' + ' + '.join(slm_name_right)
            components_equation = ' + '.join(comp_left)     + ' = ' + ' + '.join(comp_right)
            return reaction_smiles, chebi_equation, slm_id_equation, slm_name_equation, components_equation

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

        # From each right group, eliminate those participants with components that do not exist on the left side.
        right_groups_filtered: List[List[Participant]] = []
        for right_group in right_groups:
            if left_components:
                filtered_group = [participant
                                  for participant in right_group
                                  if (not participant.components)
                                  or all((comp in left_components) for comp in participant.components)]
            else:
                filtered_group = right_group

            # If one of the filtered right groups is empty, return empty list + no match.
            if not filtered_group:
                return [], 'no_component_match'
            
            right_groups_filtered.append(filtered_group)
        
        # Loop over the cartesian product of the participants in the filtered
        # right groups, skipping those with unequal components counts.
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
                # TODO: Do we keep these just for debugging? Could be expensive?
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
            iso_slm_id = rec.get('iso_slm_id')
            iso_slm_name = rec.get('iso_slm_name')
            components = tuple(rec.get('components', []) or [])
            comp_counter = rec.get('comp_counter') or Counter(components)
            comp_text = str(rec.get('Components*'))
            return Participant(chebi_id, chebi_display, side, stoich, smiles, iso_slm_id, iso_slm_name, components, comp_counter, comp_text)

        # Helper function to group participants by their class.
        def __group_participants_by_class(df_participants: pd.DataFrame) -> List[List[Participant]]:
            """
            Groups the participants (isomeric subspecies lipid ID) by their class (Rhea CHEBI ID).
            Returns List[ Rhea CHEBI ID->List[ isomeric subspecies lipid ID ] ], keeping stable order.
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

        results = []
        
        logger.info("progress_apply: Enumerating reactions")
        # Loop over the left side, grouped by Rhea ID.
        for left_reaction_side, df_left_group in tqdm(df_left.groupby('reaction_side')):
            master_id = int(left_reaction_side.split('_')[0])
            # Get the right side for this Rhea ID.
            df_right_group = df_right[df_right['MASTER_ID'] == master_id]

            signal.alarm(90000)
            try:
                # Group the hypothetical participants (isomeric subspecies lipid IDs)
                # of each side by their classes (Rhea CHEBIs).
                # At the same time, convert pandas to lightweight Python data
                # structure to make the enumeration fast.
                left_groups:  List[List[Participant]] = __group_participants_by_class(df_left_group)
                right_groups: List[List[Participant]] = __group_participants_by_class(df_right_group)

                # If any side has no groups, there is nothing to enumerate.
                if not left_groups or not right_groups:
                    continue

                # Loop over the cartesian product of the participants in the left groups:
                for left_combination in product(*left_groups):
                    # For each left combination, build the equations with the right combinations with pruning.
                    results_for_left_combination, component_match = self.__build_equations(left_combination, right_groups)
                    for reaction_smiles, chebi_equation, slm_id_equation, slm_name_equation, components_equation in results_for_left_combination:
                        if reaction_smiles:
                            results.append({
                                'MASTER_ID': master_id,
                                'chebi_equation': chebi_equation,
                                'slm_id_equation': slm_id_equation,
                                'slm_name_equation': slm_name_equation,
                                'components_equation': components_equation,
                                'reaction_smiles': reaction_smiles,
                                'component_match': component_match
                            })
            except TimeoutException:
                logger.warning("MASTER ID %s took too long and was skipped: %d",
                               master_id.split('_')[0], len(df))
                continue
            finally:
                signal.alarm(0) # Cancel alarm.

        return pd.DataFrame(results)


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
    def save_results(self, df_rhea2participants, filename, dict_rhea2bond_changes, stats):

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
        def __check_reactions_for_which_component_matching_was_impossible_with_atom_mapper(row, dict_rhea2bond_changes):
            if row['component_match'] == 'no_component_match':
                discrepancies = __summarize_discrepancies(row['reaction_smiles'])
                return len(discrepancies['bond_changes']) <= dict_rhea2bond_changes[row['MASTER_ID']]
            return True
        

        # --- Main section ---

        # self.df_rhea2participants2slm_class2iso contains only Rhea
        # participants that map to an SLM. Before we can enumerate reactions, we
        # have to add back the Rhea participants that do not map to an SLM.
        df = self.df_rhea2participants2slm_class2iso
        logger.info( "LENGTH: %8d df_rhea2participants", len(df_rhea2participants) )
        logger.info( "LENGTH: %8d df = self.df_rhea2participants2slm_class2iso", len(df) )
        logger.debug( "FORMAT: df\n%s", debug_df_first_row(df) )
        if DEBUG > 1:
            df.to_csv(os.path.join(self.output_dir, 'DEBUG_df.tsv'), sep="\t", header=True, index=False)
        
        # Get subset of df_rhea2participants columns.
        df_rhea2participants_subset = df_rhea2participants[['MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'smiles']]
        # Get subset of df_rhea2participants rows whose MASTER_ID is in df (= r2sl.df_rhea2participants2slm_class2iso).
        df_rhea2participants_subset = df_rhea2participants_subset[df_rhea2participants_subset['MASTER_ID'].isin(df['MASTER_ID'])]
        logger.info( "LENGTH: %8d df_rhea2participants_subset", len(df_rhea2participants_subset) )
        logger.debug( "FORMAT: df_rhea2participants_subset\n%s", debug_df_first_row(df_rhea2participants_subset) )

        # RIGHT merge df with df_rhea2participants_subset to add the CHEBIs from
        # df_rhea2participants_subset that have no SLM.
        logger.debug( "RIGHT merging df WITH df_rhea2participants_subset ON 'MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'smiles'" )
        df = df.merge(
            df_rhea2participants_subset, on=['MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'smiles'], how='right'
        )
        logger.info( "LENGTH: %8d df", len(df) )
        logger.debug( "FORMAT: df\n%s", debug_df_first_row(df) )
        if DEBUG > 0:
            df.to_csv(os.path.join(self.output_dir, 'DEBUG_df-after_merge.tsv'), sep="\t", header=True, index=False)
        
        df_result = self.__enumerate_reactions(df)
        stats["reactions"].append( "%8d Reactions enumerated\n" % len(df_result) )
        stats["rhea"].append( "%8d Rhea IDs (in enumerated reactions)\n" % len(df_result['MASTER_ID'].unique()) )
        if DEBUG > 1:
            logger.debug( "FORMAT: df_result\n%s", debug_df_first_row(df_result) )
            df_result.to_csv(os.path.join(self.output_dir, 'DEBUG_df_result-reactions-enumerated.tsv'), sep="\t", header=True, index=False)

        # Drop rows without reaction SMILES.
        df_result.dropna(subset = ['reaction_smiles'], inplace=True)
        stats["reactions"].append( "%8d Reactions enumerated (after filtering those without reaction SMILES)\n" % len(df_result) )
        stats["rhea"].append( "%8d Rhea IDs (after filtering those without reaction SMILES)\n" % len(df_result['MASTER_ID'].unique()) )
        
        logger.info( "progress_apply: Checking components balanced" )
        # Example input rows:
        # MASTER_ID components_equation
        # 10272 SLM:000000510 (sn1 or sn2) + nan = nan + SLM:000000510 (sn1 or sn2) + nan
        # 11488	nan + SLM:000000510 (free fatty acid) = SLM:000000510 (acyl) + nan
        parsed = df_result['components_equation'].progress_apply(__get_component_changes).apply(pd.Series)
        df_result = pd.concat([df_result, parsed], axis=1)
        if DEBUG > 0:
            df_result.to_csv(os.path.join(self.output_dir, 'DEBUG_df_result-before-applying-components-balanced.tsv'), sep="\t", header=True, index=False)
        # Store exceptions:
        # 63596, 78043, 77759, 78435 have triglyceride twice in reactants, therefore sn1 is not unique, 
        # there is sn1 of the first tryglyceride, and sn2 of the second, which makes the code confused
        df_change_of_sn_by_design = df_result[df_result['MASTER_ID'].isin([63596, 78043, 77759, 78435])]
        df_change_of_sn_by_design['component_match'] = 'no_component_match'
        # Get list of boolean '_change' columns and remove rows where any of them is True.
        change_cols = [col for col in df_result.columns if col.endswith('_change')]
        df_result = df_result[~df_result[change_cols].any(axis=1)]
        stats["reactions"].append( "%8d Reactions enumerated (after filtering unbalanced components)\n" % len(df_result) )
        stats["rhea"].append( "%8d Rhea IDs (after filtering unbalanced components)\n" % len(df_result['MASTER_ID'].unique()) )
        # Now add back the exceptions.
        df_result = pd.concat([df_result, df_change_of_sn_by_design])
        if DEBUG > 0:
            df_result.to_csv(os.path.join(self.output_dir, 'DEBUG_df_result-after-checking-components-balanced.tsv'), sep="\t", header=True, index=False)
        stats["reactions"].append( "%8d Reactions enumerated (after adding back triglyceride exceptions)\n" % len(df_result) )
        stats["rhea"].append( "%8d Rhea IDs (after adding back triglyceride exceptions)\n" % len(df_result['MASTER_ID'].unique()) )

        # Filter unbalanced reactions.
        rxn = Reaction()
        logger.info( "progress_apply: Checking reaction balanced" )
        df_result['balanced'] = df_result['reaction_smiles'].progress_apply(
            lambda reaction_smiles: rxn.check_reaction_balance(reaction_string=reaction_smiles, format="smiles", account_H=True)
        )
        df_result = df_result[df_result['balanced'] == True]
        df_result.drop(columns=['balanced'], inplace=True)
        stats["reactions"].append( "%8d Reactions enumerated (after filtering unbalanced reactions)\n" % len(df_result) )
        stats["rhea"].append( "%8d Rhea IDs (after filtering unbalanced reactions)\n" % len(df_result['MASTER_ID'].unique()) )
        
        # Filter reactions with too many bond changes.
        logger.info( "progress_apply: Checking number of bond changes" )
        df_result['bond_breakage_max_4'] = df_result.progress_apply(__check_reactions_for_which_component_matching_was_impossible_with_atom_mapper, axis=1, args=[dict_rhea2bond_changes,])
        df_result = df_result[df_result['bond_breakage_max_4'] == True]
        stats["reactions"].append( "%8d Reactions enumerated (after filtering too many bond changes)\n" % len(df_result) )
        stats["rhea"].append( "%8d Rhea IDs (after filtering too many bond changes)\n" % len(df_result['MASTER_ID'].unique()) )

        # Generate RInChI and Web-RInChIKey.
        logger.info( "progress_apply: Generating RInChI" )
        rinchi = RInChI()
        if len(df_result) > 0:
            df_result[['RInChI', 'Web-RInChIKey']] = df_result.progress_apply(
                lambda x: rinchi.error_handle_wrap_rinchi(x['reaction_smiles']),
                axis=1, result_type='expand'
            )
        
            # Filter reactions without Web-RInChIKey (without structures) and duplicate Web-RInChIKey.
            df_result.dropna(subset=['Web-RInChIKey'], inplace=True)
            df_result.drop_duplicates(subset=['MASTER_ID', 'Web-RInChIKey', 'reaction_smiles'], inplace=True)
            stats["reactions"].append( "%8d Reactions enumerated (after filtering empty and duplicated Web-RInChIKey)\n" % len(df_result) )
            stats["rhea"].append( "%8d Rhea IDs (after filtering empty and duplicated Web-RInChIKey)\n" % len(df_result['MASTER_ID'].unique()) )

            # Save final result.
            df_result.to_csv(os.path.join(self.output_dir, filename), sep='\t', index=False)
