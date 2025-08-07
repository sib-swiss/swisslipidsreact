import os
import re
import pandas as pd
import numpy as np
from itertools import product
import signal
from rxnmapper import RXNMapper
from rdkit import Chem

from pyrheadb.Reaction import Reaction
from pyrheadb.RInChI import RInChI

from tqdm import tqdm

from .FA_lists import positions
tqdm.pandas()

class TimeoutException(Exception):
    pass

def handler(signum, frame):
    raise TimeoutException()

# Set the signal handler
signal.signal(signal.SIGALRM, handler)


class RheaToSwisslipidsDf():
    
    def __init__(self, output_dir='', timestamp='now'):
        self.output_dir=output_dir
        self.timestamp=timestamp

    def get_df_swiss_lipids_chebi_rhea(self, df_swiss_lipids_chebi, df_rhea_chebi):
        self.df_swiss_lipids_chebi_rhea = df_rhea_chebi.merge(df_swiss_lipids_chebi, on='chebi_id', how='inner')

        # ---------- Merge to Obtain Rhea Reactions of Specific Subspecies ----------

    def get_df_rhea_descendant(self, rhea_lipid_to_descendant_df, df_swisslipids):
        df_rhea_descendant = self.df_swiss_lipids_chebi_rhea.merge(
            rhea_lipid_to_descendant_df, left_on='Lipid ID', right_on='rhea_lipid_id', how='left'
        )

        df_rhea_descendant = df_rhea_descendant[[
            'MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef',
            'rhea_lipid_id', 'isomeric_subspecies_descendant_lipid_id', 'smiles'
        ]]

        df_rhea_descendant = df_rhea_descendant.merge(
            df_swisslipids, left_on='isomeric_subspecies_descendant_lipid_id',
            right_on='Lipid ID', how='left'
        )

        df_rhea_descendant.drop(columns=['Lipid ID', 'Level', 'Lipid class*'], inplace=True)
        # print('Rhea * swisslipid isomeric subspecies', len(df_rhea_descendant))
        # print('Unique swisslipid isomeric subspecies in Rhea', len(set(df_rhea_descendant['isomeric_subspecies_descendant_lipid_id'])))

        self.df_rhea_descendant = df_rhea_descendant
        return len(df_rhea_descendant), len(set(df_rhea_descendant['isomeric_subspecies_descendant_lipid_id']))


    def extract_stoichiometric_components(self, row):

        # Only include columns if they exist in the DataFrame
        columns_to_extract = positions
        components = []

        stoich = int(row['stoich_coef']) if pd.notnull(row['stoich_coef']) else 1  # Default to 1 if missing

        for col in columns_to_extract:
            if col in row and pd.notnull(row[col]) and isinstance(row[col], str) and row[col].strip():
                components.extend([row[col]] * stoich)

        return components
    
    # ---------- Reaction Generation Functions ----------
    def build_equations(self, combination, df_right):

        res_final = []
        res_backup = []
        """Construct reaction SMILES and corresponding equations from a reaction component combination."""
        df_left = pd.DataFrame(combination)

        def get_smiles(row):
            """Select SMILES (pH7.3) if available, otherwise use fallback SMILES."""
            return row['SMILES (pH7.3)'] if pd.notna(row['SMILES (pH7.3)']) else row['smiles']
        
        def combination_to_equations(df_left, df_right_temp):
            combination = pd.concat([df_left, df_right_temp], ignore_index=True)

            left, right = [], []
            chebi_eq_parts_left, chebi_eq_parts_right = [], []
            swisslipids_eq_parts_left, swisslipids_eq_parts_right = [], []
            components_parts_left, components_parts_right = [], []

            for index, row in combination.iterrows():
                smiles = get_smiles(row)
                repeated_smiles = [smiles] * int(row['stoich_coef'])

                # Equations based on ChEBI or SwissLipids
                chebi_term = f"{int(row['stoich_coef'])} {row['CHEBI']}" if pd.notna(row['CHEBI']) else f"{int(row['stoich_coef'])} {row['chebi_id']}"
                swiss_term = f"{int(row['stoich_coef'])} {row['isomeric_subspecies_descendant_lipid_id']}" if pd.notna(row['isomeric_subspecies_descendant_lipid_id']) else f"{int(row['stoich_coef'])} NA"

                components = str(row['Components*'])

                if row['reaction_side'].endswith('L'):
                    left.extend(repeated_smiles)
                    chebi_eq_parts_left.append(chebi_term)
                    swisslipids_eq_parts_left.append(swiss_term)
                    components_parts_left.append(components)
                else:
                    right.extend(repeated_smiles)
                    chebi_eq_parts_right.append(chebi_term)
                    swisslipids_eq_parts_right.append(swiss_term)
                    components_parts_right.append(components)
            
            reaction = '.'.join(left) + '>>' + '.'.join(right)
            chebi_equation = ' + '.join(chebi_eq_parts_left) + ' = ' + ' + '.join(chebi_eq_parts_right)
            swisslipids_equation = ' + '.join(swisslipids_eq_parts_left) + ' = ' + ' + '.join(swisslipids_eq_parts_right)
            components_equation = ' + '.join(components_parts_left) + ' = ' + ' + '.join(components_parts_right)

            return reaction, chebi_equation, swisslipids_equation, components_equation

        components_left = [
                x
                for xs in df_left['components']
                for x in xs
            ]

        if components_left:
            for sn in positions:
                df_right = df_right[df_right[sn].isin(components_left+[np.nan, ''])]
        
        components_right = [
                x
                for xs in df_right['components']
                for x in xs
            ]

        chebi_groups = [grp.to_dict('records') for _, grp in df_right.groupby('chebi_id')]
        
        if len(df_right)==0:
            #print(df_left)
            return [], 'no_component_match'
        for combo in product(*chebi_groups):
            df_right_temp = pd.DataFrame(combo)

            components_right = [
                x
                for xs in df_right_temp['components']
                for x in xs
            ]

            components_left.sort()
            components_right.sort()
            if ','.join(components_left) != ','.join(components_right) or len(components_left)==0 and len(components_right)==0:
                res_final.append([None, None, None, None])
                res_backup.append(combination_to_equations(df_left, df_right_temp))
            
            else:
                res_final.append(combination_to_equations(df_left, df_right_temp))
        res_final = [li for li in res_final if li[0] is not None]

        if res_final:
            return res_final, 'yes_component_match'
        return res_backup, 'no_component_match'

    def generate_combinations_and_reactions(self, df):
        """Generates all valid reaction SMILES and equation strings for each MASTER_ID."""

        print('Extracting stoichiometric coefficients per compound ->')
        df['components'] = df.progress_apply(self.extract_stoichiometric_components, axis=1)
        all_results = []
        df['side'] = df['reaction_side'].apply(lambda x: x.split('_')[1])
        df_left = df[df['side']=='L']
        df_right = df[df['side']=='R']

        total_gen_attempted = len(df_left['reaction_side'].unique())

        print('Generating combinations of reactants and products ->')
        for reaction_side_master_id, group in tqdm(df_left.groupby('reaction_side')):
            df_right_master_id = df_right[df_right['MASTER_ID']==int(reaction_side_master_id.split('_')[0])]
            signal.alarm(90000)
            try:
                chebi_groups = [grp.to_dict('records') for _, grp in group.groupby('chebi_id')]
                for combo in product(*chebi_groups):
                    res_total, component_match = self.build_equations(combo, df_right_master_id)
                    for reaction, chebi_eq, swisslipids_eq, components_equation in res_total:
                        if reaction:
                            all_results.append({
                                'MASTER_ID': int(reaction_side_master_id.split('_')[0]),
                                'chebi_equation': chebi_eq,
                                'swisslipids_equation': swisslipids_eq, 
                                'components_equation': components_equation,
                                'reaction_smiles': reaction,
                                'component_match':component_match
                            })
            except TimeoutException:
                print(f"MASTER ID {reaction_side_master_id.split('_')[0]} took too long and was skipped.", len(df))
                continue
            finally:
                signal.alarm(0)  # Cancel alarm
        return pd.DataFrame(all_results), total_gen_attempted

    def extract_atom_environment(self, mol):
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

    def summarize_discrepancies(self, reaction_smiles):
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

        reactant_smiles, product_smiles = mapped_rxn.split(">>")
        reactants = [Chem.MolFromSmiles(s) for s in reactant_smiles.split(".")]
        products  = [Chem.MolFromSmiles(s) for s in product_smiles.split(".")]

        react_env = {}
        for mol in reactants:
            react_env.update(self.extract_atom_environment(mol))

        prod_env = {}
        for mol in products:
            prod_env.update(self.extract_atom_environment(mol))

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

    def check_reactions_for_which_component_matching_was_impossible_with_atom_mapper(self, row):
        if row['component_match']=='no_component_match':
            discrepancies = self.summarize_discrepancies(row['reaction_smiles'])
            return len(discrepancies['bond_changes'])<5
        return True

    def save_results(self, rheadf, df, rhea_reactions, filename):

        df_rhea_cut = rheadf[['MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'smiles']]
        df_rhea_cut = df_rhea_cut[df_rhea_cut['MASTER_ID'].isin(df['MASTER_ID'])]
        df = df.merge(
            df_rhea_cut, on=['MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'smiles'], how='right'
        )
        result_df, total_gen_attempted = self.generate_combinations_and_reactions(df)
        rhea_reactions['5_after_enumerating_combinations'] = rhea_reactions['MASTER_ID'].isin(result_df['MASTER_ID'])
        result_df.dropna(subset = ['reaction_smiles'], inplace=True)
        rhea_reactions['6_after_dropping_na_reaction_smiles'] = rhea_reactions['MASTER_ID'].isin(result_df['MASTER_ID'])
        MASTER_ID_for_specific_fatty_acids_before_filtering_out_the_unbalanced = len(set(result_df['MASTER_ID']))
        num_reactions_to_check_for_balance = len(result_df) # Total # of reactions to check for balance
        
        def extract_position_map(eq_side: str) -> dict[str, str]:
            """
            Given one side of the component equation (left or right), return {position: SLM}.
            Uses the same parsing logic from the extractall + explode + pivot code.
            """
            # pattern like "SLM:000000123 (sn1 or sn2)"
            pattern = r'(SLM:\d+)\s+\(([^)]+)\)'
            matches = re.findall(pattern, eq_side)

            extracted_rows = []
            for slm, pos_string in matches:
                positions_split = re.split(r'\s*(?:or|and|,)\s*', pos_string)
                for pos in positions_split:
                    pos = pos.strip()
                    if pos in positions:
                        extracted_rows.append((pos, slm))

            # Keep only one SLM per position
            pos_map = {}
            for pos, slm in extracted_rows:
                if pos not in pos_map:  # take first one only
                    pos_map[pos] = slm
            return pos_map

        def summarise(eq: str):
            """For a given equation string, return a dict with position summaries and change flags."""
            left, right = eq.split('=', 1)
            left_map = extract_position_map(left)
            right_map = extract_position_map(right)

            out = {}
            for p in positions:
                start = left_map.get(p)
                end = right_map.get(p)

                if start is None and end is None:
                    out[p] = ''
                    out[f'{p}_change'] = False
                else:
                    out[p] = f'{start or "None"}-{end or "None"}'
                    out[f'{p}_change'] = (start is not None and end is not None and start != end)
            return out

        print('Checking component balance')
        parsed = result_df['components_equation'].progress_apply(summarise).apply(pd.Series)
        result_df = pd.concat([result_df, parsed], axis=1)

        # List all columns ending with "_change"
        change_cols = [col for col in result_df.columns if col.endswith('_change')]

        # Select rows where at least one position has changed
        # df_changed = result_df[result_df[change_cols].any(axis=1)]
        # df_changed.drop_duplicates(subset=['MASTER_ID'],inplace=True)
        # df_changed.to_csv('changes_examples.tsv', sep='\t', index=False)

        result_df = result_df[~result_df[change_cols].any(axis=1)]

        rhea_reactions['7_after_dropping_incorrect_component_matches'] = rhea_reactions['MASTER_ID'].isin(result_df['MASTER_ID'])

        # Filter for balanced reactions
        rxn = Reaction()
        print('checking reaction balance')
        result_df['balanced'] = result_df['reaction_smiles'].progress_apply(rxn.check_reaction_balance)
        result_df = result_df[result_df['balanced']==True]
        result_df.drop(columns=['balanced'], inplace=True)
        MASTER_ID_for_specific_fatty_acids_after_filtering_out_the_unbalanced = len(set(result_df['MASTER_ID']))
        rhea_reactions['8_after_dropping_unbalanced_reactions'] = rhea_reactions['MASTER_ID'].isin(result_df['MASTER_ID'])
        num_reactions_after_filtering_out_the_unbalanced = len(result_df)
        self.rxn_mapper = RXNMapper()
        print('Check that less than 5 bonds are broken at the same time')
        result_df['bond_breakage_less_than_5'] = result_df.progress_apply(self.check_reactions_for_which_component_matching_was_impossible_with_atom_mapper, axis=1)
        result_df = result_df[result_df['bond_breakage_less_than_5']==True]
        rhea_reactions.to_csv(os.path.join(self.output_dir, f'{self.timestamp}_rhea_reactions_overview.tsv'), sep='\t', index=False)

        # Generate RInChI and Web-RInChIKey
        rinchi = RInChI()
        print("Generating RInChI")
        if len(result_df)>0:
            result_df[['RInChI', 'Web-RInChIKey']] = result_df.progress_apply(
                lambda x: rinchi.error_handle_wrap_rinchi(x['reaction_smiles']),
                axis=1, result_type='expand'
            )
        
        # Clean the data of the reactions withour Web-RInChIKey (without structures) and dulpicates
        result_df.dropna(subset=['Web-RInChIKey'], inplace=True)
        result_df.drop_duplicates(subset=['MASTER_ID', 'Web-RInChIKey', 'reaction_smiles'], inplace=True)
        
        # Save final result
        result_df.to_csv(os.path.join(self.output_dir, filename), sep='\t', index=False)

        return total_gen_attempted, \
              MASTER_ID_for_specific_fatty_acids_before_filtering_out_the_unbalanced, \
               num_reactions_to_check_for_balance, \
               MASTER_ID_for_specific_fatty_acids_after_filtering_out_the_unbalanced, \
               num_reactions_after_filtering_out_the_unbalanced
