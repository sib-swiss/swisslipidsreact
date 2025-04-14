import os
import pandas as pd
from itertools import product

from pyrheadb.Reaction import Reaction
from pyrheadb.RInChI import RInChI

class RheaToSwisslipidsDf():
    
    def __init__(self):
        pass

    def get_df_swiss_lipids_chebi_rhea(self, df_swiss_lipids_chebi, df_rhea_chebi):
        self.df_swiss_lipids_chebi_rhea = df_swiss_lipids_chebi.merge(df_rhea_chebi, on='chebi_id', how='inner')
        print('# swiss lipids * chebi * rhea', len(self.df_swiss_lipids_chebi_rhea))

        # ---------- Merge to Obtain Rhea Reactions of Specific Subspecies ----------

    def get_df_rhea_descendant(self, rhea_lipid_to_descendant_df, df_swisslipids):
        df_rhea_descendant = self.df_swiss_lipids_chebi_rhea.merge(
            rhea_lipid_to_descendant_df, left_on='Lipid ID', right_on='rhea_lipid_id', how='inner'
        )

        df_rhea_descendant = df_rhea_descendant[[
            'MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef',
            'rhea_lipid_id', 'isomeric_subspecies_descendant_lipid_id'
        ]]

        df_rhea_descendant = df_rhea_descendant.merge(
            df_swisslipids, left_on='isomeric_subspecies_descendant_lipid_id',
            right_on='Lipid ID', how='inner'
        )
        df_rhea_descendant.drop(columns=['Lipid ID', 'Level', 'Lipid class*'], inplace=True)
        print('Rhea * swisslipid isomeric subspecies', len(df_rhea_descendant))
        print('Unique swisslipid isomeric subspecies in Rhea', len(set(df_rhea_descendant['isomeric_subspecies_descendant_lipid_id'])))

        self.df = df_rhea_descendant

    def filter_only_specific_FA_set(self, rheadf, FA_SLM_list=[]):
        """
        Keep only the part of the dataframe related to particular FAs
        """

        filtered_df_rhea_descendant = self.df.copy()
        # filtered_df_rhea_descendant.dropna(subset=['Components*'], inplace=True)
        filtered_df_rhea_descendant['Components*'] = filtered_df_rhea_descendant['Components*'].astype(str)

        # Only retain rows where only defined fatty acids are present
        filtered_df_rhea_descendant['only_defined_FA'] = filtered_df_rhea_descendant['Components*'].apply(
            self.if_only_defined_FA, FA_SLMs=FA_SLM_list
        )
        filtered_df_rhea_descendant = filtered_df_rhea_descendant[
            filtered_df_rhea_descendant['only_defined_FA']
        ]
        print("len(filtered_df_rhea_descendant)", len(filtered_df_rhea_descendant))

        df_rhea_cut = rheadf[['MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef', 'smiles']]
        df_rhea_cut = df_rhea_cut[df_rhea_cut['MASTER_ID'].isin(filtered_df_rhea_descendant['MASTER_ID'])]

        filtered_df_rhea_descendant = filtered_df_rhea_descendant.merge(
            df_rhea_cut, on=['MASTER_ID', 'reaction_side', 'chebi_id', 'stoich_coef'], how='right'
        )
        filtered_df_rhea_descendant.to_csv(os.path.join('..','..', 'results','filtered_df_rhea_descendant.tsv'), sep='\t', index=False)
        self.filtered_df_rhea_descendant = filtered_df_rhea_descendant

    def if_only_defined_FA(self, components_string, FA_SLMs):
        """
        Returns True if the components string contains only the defined fatty acids (SLMs).
        Used to filter for lipids that only contain specific fatty acids (e.g., palmitic acid).
        """
        for SLM in FA_SLMs:
            components_string = components_string.replace(SLM, '')
        return 'SLM' not in components_string
    
    # ---------- Reaction Generation Functions ----------
    def build_equations(self, combination):
        """Construct reaction SMILES and corresponding equations from a reaction component combination."""

        def get_smiles(row):
            """Select SMILES (pH7.3) if available, otherwise use fallback SMILES."""
            return row['SMILES (pH7.3)'] if pd.notna(row['SMILES (pH7.3)']) else row['smiles']
    
        left, right = [], []
        chebi_eq_parts_left, chebi_eq_parts_right = [], []
        swisslipids_eq_parts_left, swisslipids_eq_parts_right = [], []
        components_parts_left, components_parts_right = [], []

        for row in combination:
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

    def generate_combinations_and_reactions(self, df):
        """Generates all valid reaction SMILES and equation strings for each MASTER_ID."""
        all_results = []
        for master_id, group in df.groupby('MASTER_ID'):
            chebi_groups = [grp.to_dict('records') for _, grp in group.groupby('chebi_id')]
            for combo in product(*chebi_groups):
                reaction, chebi_eq, swisslipids_eq, components_equation = self.build_equations(combo)
                all_results.append({
                    'MASTER_ID': master_id,
                    'chebi_equation': chebi_eq,
                    'swisslipids_equation': swisslipids_eq, 
                    'components_equation': components_equation,
                    'reaction_smiles': reaction
                })
        return pd.DataFrame(all_results)

    # ---------- Generate, Filter, and Save Final Results ----------
    def save_results_per_fa(self,filename):
        self.save_results(self.filtered_df_rhea_descendant, filename)

    def save_results(self, df, filename):
        result_df = self.generate_combinations_and_reactions(df)
        print('# MASTER ID for specific fatty acids before filtering out the unbalanced', len(set(result_df['MASTER_ID'])))

        # Filter for balanced reactions
        rxn = Reaction()
        result_df['balanced'] = result_df['reaction_smiles'].apply(rxn.check_reaction_balance)
        result_df = result_df[result_df['balanced'] == True]
        print('# MASTER ID for specific fatty acids after filtering out the unbalanced', len(set(result_df['MASTER_ID'])))

        # Generate RInChI and Web-RInChIKey
        rinchi = RInChI()
        result_df[['RInChI', 'Web-RInChIKey']] = result_df.apply(
            lambda x: rinchi.error_handle_wrap_rinchi(x['reaction_smiles']),
            axis=1, result_type='expand'
        )

        # Save final result
        result_df.drop(columns=['balanced'], inplace=True)
        result_df.to_csv(os.path.join('..', '..', 'results', filename), sep='\t', index=False)
