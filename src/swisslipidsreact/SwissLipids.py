import pandas as pd
import os
import networkx as nx
import importlib.resources
from platformdirs import user_cache_dir

from .FA_lists import positions, get_FA_list, FAS_15, FAS_85, FAS_79, FOH_15, PAL_C16, PALOH_C16, PAL_C16_OCT_C18, SPHINGO_23

flag_fast_exec = True

class SwissLipids():

    def __init__(self, output_dir='', timestamp='now'):
        # get a safe cache directory
        self.cache_dir = user_cache_dir("swisslipidsreact")
        os.makedirs(self.cache_dir, exist_ok=True)
        self.output_dir = output_dir
        self.timestamp=timestamp

    def preporcess_swisslipids(self):
        """
        This function adds (free fatty acid) as component to FA and FA-COA and n-acyls to acylethanolamines
        """

        with importlib.resources.files("swisslipidsreact.package_data").joinpath("lipids.tsv").open("rb") as f:
            self.swisslipids = pd.read_csv(f, sep="\t", encoding="latin-1",
                usecols=['Name', 'Lipid ID', 'CHEBI', 'Level', 'Lipid class*', 'Components*', 'SMILES (pH7.3)'],
                dtype={'Lipid ID': str, 'CHEBI': str, 'Level': str, 'Lipid class*': str,
                    'Components*': str, 'SMILES (pH7.3)': str})

        # assign level to be isomeric subspecies for all swisslipids that are not class lipids (compounds with *)
        self.swisslipids.loc[~self.swisslipids['SMILES (pH7.3)'].isna() & 
        self.swisslipids['Level'].isna() & ~self.swisslipids['SMILES (pH7.3)'].str.contains('*', regex=False, na=False),
        'Level'] = 'Isomeric subspecies'

        # Normalize prime characters to ASCII apostrophe (')
        prime_variants = ['′', 'ʹ', '´']  # U+2032, U+02B9, U+00B4
        for variant in prime_variants:
            self.swisslipids['Components*'] = self.swisslipids['Components*'].str.replace(variant, "'", regex=False)

        # Function to clean Components* column
        def clean_components(row):
            """
            remove class from the component
            """
            components = row['Components*']
            lipid_class = row['Lipid class*']
            
            if pd.isna(components):
                return components  # Leave NaN as is

            # Split and filter
            filtered = [comp.strip() for comp in components.split(' / ') if lipid_class not in comp]
            return ' / '.join(filtered) if filtered else None

        # Apply clean_components function
        self.swisslipids['Components*'] = self.swisslipids.apply(clean_components, axis=1)

        # Filter to retain only isomeric subspecies in a separate df
        self.get_isomeric_subspecies_table() # generates self.df_isomeric_subspecies
        
        self.get_lipid_class_graph()

        # Assign free fatty acid as component for all FA and FA-CoA
        FA = "SLM:000000984"
        F_alcohol = "SLM:000390053"

        FAs = self.get_lipid_to_descendant_df([FA])['Lipid ID'].to_list()
        F_alcohols = self.get_lipid_to_descendant_df([F_alcohol])['Lipid ID'].to_list()
        
        # Assign (free fatty acid) to Components* for FAs
        self.swisslipids.loc[self.swisslipids['Lipid ID'].isin(FAs), 'Components*'] = (
        self.swisslipids['Lipid ID'] + ' (free fatty acid)' )

        # Assign (free fatty alcohol) to Components* for FAs
        self.swisslipids.loc[self.swisslipids['Lipid ID'].isin(F_alcohols), 'Components*'] = (
        self.swisslipids['Lipid ID'] + ' (free fatty alcohol)' )

        # helper to avoid repeating:
        def load_package_tsv(module, filename, sep="\t"):
            with importlib.resources.files(module).joinpath(filename).open("r") as f:
                return pd.read_csv(f, sep=sep)

        # Now load the files:
        df_facoa_to_ffa = load_package_tsv(
            "swisslipidsreact.package_data.components_correction", 
            "free fatty acids per compound.tsv"
        )
        comp_to_ffa_dict = dict(zip(df_facoa_to_ffa['Lipid ID'], df_facoa_to_ffa['free fatty acid']))

        df_comp_to_nacyl_dict = load_package_tsv(
            "swisslipidsreact.package_data.components_correction",
            "nacyls per compound.tsv"
        )
        comp_to_nacyl_dict = dict(zip(df_comp_to_nacyl_dict['Lipid ID'], df_comp_to_nacyl_dict['n-acyl']))

        df_comp_to_sn1_dict = load_package_tsv(
            "swisslipidsreact.package_data.components_correction",
            "sn1 per compound.tsv"
        )
        comp_to_sn1_dict = dict(zip(df_comp_to_sn1_dict['Lipid ID'], df_comp_to_sn1_dict['sn1']))

        df_comp_to_sn2_dict = load_package_tsv(
            "swisslipidsreact.package_data.components_correction",
            "sn2 per compound.tsv"
        )
        comp_to_sn2_dict = dict(zip(df_comp_to_sn2_dict['Lipid ID'], df_comp_to_sn2_dict['sn2']))
        
        def assign_components(row, comp_to_ffa_dict, comp_to_nacyl_dict, comp_to_sn1_dict, comp_to_sn2_dict):
            lipid_id = row['Lipid ID']
            nacyl_part = None
            acid_part = None
            sn1_part = None
            sn2_part = None
            # Each part can only be mentioned once in the components string, so I make logic based on this
            if lipid_id in comp_to_nacyl_dict:
                nacyl_part = f"{comp_to_nacyl_dict[lipid_id]} (n-acyl)"
            if lipid_id in comp_to_ffa_dict:
                acid_part = f"{comp_to_ffa_dict[lipid_id]} (free fatty acid)"
            if lipid_id in comp_to_sn1_dict:
                sn1_part = f"{comp_to_sn1_dict[lipid_id]} (sn1)"
            if lipid_id in comp_to_sn2_dict:
                sn2_part = f"{comp_to_sn2_dict[lipid_id]} (sn2)"
            
            if sn1_part and sn2_part:
                return sn1_part + ' / ' + sn2_part
            elif nacyl_part:
                if type(row['Components*'])==str:
                    return row['Components*'] + nacyl_part
                return nacyl_part
            elif acid_part:
                return acid_part
            return row['Components*']

        # Assign values to Components* for CoAs
        self.swisslipids['Components*'] = self.swisslipids.apply(assign_components, axis=1, args = [comp_to_ffa_dict, comp_to_nacyl_dict, comp_to_sn1_dict, comp_to_sn2_dict, ])

        """
                Class	Component 	Component 
                            N-acyl	O-acyl
		SLM:000389330	2-[(9Z)-octadecenoylamino]ethyl (9Z)-octadecenoate	SLM:000508875	SLM:000000418	SLM:000000418
        SLM:000389331	2-[(5Z,8Z,11Z,14Z)-eicosatetraenoylamino]ethyl (9Z)-octadecenoate	SLM:000508875	SLM:000000296	SLM:000000418
        """
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389330', 'Components*'] = 'SLM:000000418 (n-acyl) / SLM:000000418 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389331', 'Components*'] = 'SLM:000000296 (n-acyl) / SLM:000000418 (acyl)'

        """
            		Class	Component 	Component 
			            acyl	N-acyl
        SLM:000389332	1-(9Z)-octadecenoyl-N-(acetyl)-sphing-4-enine	SLM:000501524	SLM:000000418	SLM:000000449
        SLM:000389333	1-hexadecanoyl-N-(acetyl)-sphing-4-enine	SLM:000501524	SLM:000000510	SLM:000000449
        SLM:000389335	1-(5Z,8Z,11Z,14Z)-eicosatetraenoyl-N-(acetyl)-sphing-4-enine	SLM:000501524	SLM:000000296	SLM:000000449
        SLM:000389825	1-octadecanoyl-N-(acetyl)-sphing-4-enine	SLM:000501524	SLM:000000826	SLM:000000449
        SLM:000389826	1-(9Z,12Z)-octadecadienoyl-N-(acetyl)-sphing-4-enine	SLM:000501524	SLM:000000407	SLM:000000449
        SLM:000389982	1-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-N-(acetyl)-sphing-4-enine	SLM:000501524	SLM:000001084	SLM:000000449
        """
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389332', 'Components*'] = 'SLM:000000418 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389333', 'Components*'] = 'SLM:000000510 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389335', 'Components*'] = 'SLM:000501524 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389825', 'Components*'] = 'SLM:000000826 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389826', 'Components*'] = 'SLM:000000407 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389982', 'Components*'] = 'SLM:000001084 (acyl)'

        """
                        Class	Component 
                    acyl
        SLM:000000509	All-trans-retinyl hexadecanoate	SLM:000000982	SLM:000000510
        SLM:000389419	All-trans-retinyl (9Z)-octadecenoate	SLM:000000982	SLM:000000418
        SLM:000389822	All-trans-retinyl octadecanoate	SLM:000000982	SLM:000000826
        SLM:000598073	all-trans-retinyl heptanoate	SLM:000000982	SLM:000389947
        SLM:000598105	all-trans-retinyl octanoate	SLM:000000982	SLM:000000853
        SLM:000598106	all-trans-retinyl dodecanoate	SLM:000000982	SLM:000000719
        """
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000000509', 'Components*'] = 'SLM:000000510 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389419', 'Components*'] = 'SLM:000000418 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000389822', 'Components*'] = 'SLM:000000826 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000598073', 'Components*'] = 'SLM:000389947 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000598105', 'Components*'] = 'SLM:000000853 (acyl)'
        self.swisslipids.loc[self.swisslipids['Lipid ID'] == 'SLM:000598106', 'Components*'] = 'SLM:000000719 (acyl)'

        # Add extra rows with CoAs

        with importlib.resources.files("swisslipidsreact.package_data").joinpath("lipids_7rows_extra_CoAs.tsv").open("r") as f:
            df_extra_rows_coas = pd.read_csv(f, sep="\t")

        self.swisslipids = pd.concat([self.swisslipids, df_extra_rows_coas], ignore_index=True)

        # your cache file path
        cache_file = os.path.join(self.cache_dir, "lipids_preprocessed.tsv")
        self.swisslipids.to_csv(cache_file, sep="\t", index=False)

    # ---------- Load SwissLipids Data ----------
    def read_swisslipids_from_file(self):

        if not os.path.exists(os.path.join(self.cache_dir, "lipids_preprocessed.tsv")):
            self.preporcess_swisslipids()

        lipids_components_split_cache_file = os.path.join(self.cache_dir, 'lipids_components_split.tsv')
        if not os.path.exists(lipids_components_split_cache_file) or flag_fast_exec==False:
            self.swisslipids = pd.read_csv(
                os.path.join(self.cache_dir, "lipids_preprocessed.tsv"), sep='\t', encoding='latin-1',
                usecols=['Name', 'Lipid ID', 'CHEBI', 'Level', 'Lipid class*', 'Components*', 'SMILES (pH7.3)'],
                dtype={'Lipid ID': str, 'CHEBI': str, 'Level': str, 'Lipid class*': str,
                    'Components*': str, 'SMILES (pH7.3)': str}
            )

            self.swisslipids.loc[~self.swisslipids['SMILES (pH7.3)'].isna() & 
            self.swisslipids['Level'].isna() & ~self.swisslipids['SMILES (pH7.3)'].str.contains('*', regex=False, na=False),
            'Level'] = 'Isomeric subspecies'

            print('splitting the components into the positions')
            # Define possible positions

            # Initialize columns to None
            for pos in positions:
                self.swisslipids[pos] = None

            # Work only on rows where Level == 'Isomeric subspecies'
            mask = self.swisslipids['Level'] == 'Isomeric subspecies'
            df_iso = self.swisslipids[mask].copy()

            # Extract (SLM ID, position info)
            pattern = r'(SLM:\d+)\s+\(([^)]+)\)'
            extracted = df_iso['Components*'].str.extractall(pattern)
            extracted.columns = ['SLM', 'positions']
            extracted['row'] = extracted.index.get_level_values(0)

            # Expand position info (e.g. "sn1 or sn2")
            extracted = extracted.assign(position=extracted['positions'].str.split(r'\s*(?:or|,|and)\s*')).explode('position')
            extracted['position'] = extracted['position'].str.strip()

            # Keep only defined positions
            extracted = extracted[extracted['position'].isin(positions)]
            if len(extracted[~extracted['position'].isin(positions)]>0):
                print(extracted[~extracted['position'].isin(positions)])

            # Drop duplicates to enforce single SLM per position
            extracted = extracted.drop_duplicates(subset=['row', 'position'])

            # Pivot to wide format
            components_wide = extracted.pivot(index='row', columns='position', values='SLM')

            # Assign values back to original df only for Isomeric subspecies rows
            for pos in positions:
                self.swisslipids.loc[mask, pos] = df_iso.index.map(components_wide.get(pos))
            self.swisslipids.to_csv(lipids_components_split_cache_file, sep='\t', index=False)
            print('DONE: splitting the components into the positions')
        else:
            self.swisslipids = pd.read_csv(lipids_components_split_cache_file, sep='\t', low_memory=False)
        
        # Filter to retain only isomeric subspecies in a seeparate df
        self.get_isomeric_subspecies_table() # generates self.df_isomeric_subspecies
            
        self.get_lipid_class_graph()
    
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
    

    def filter_curated_biologically_relevant_isomeric_subspecies_only(self, curated_fa_list_run=True):
        """
        Creates df prefiltered to only have lipids specific for a particular positions / FA
        Columns:
        ['rhea_lipid_id', 'isomeric_subspecies_descendant_lipid_id', 'Lipid ID',
       'Level', 'Name', 'Lipid class*', 'Components*', 'SMILES (pH7.3)',
       'CHEBI', 'sn1'', 'sn2'', 'n-acyl', 'sn1', 'sn2', 'sn3',
       'free fatty acid', 'free fatty alcohol']
        """

        # Generate FA pools only once
        pos_descr_to_FA_list = {
        'FOH_15': self.SLMs_from_CHEBIs(get_FA_list(FOH_15)),
        'SPHINGO_23': self.SLMs_from_CHEBIs(get_FA_list(SPHINGO_23)),
        'PAL_C16': self.SLMs_from_CHEBIs(get_FA_list(PAL_C16)),
        'PAL_C16_OCT_C18': self.SLMs_from_CHEBIs(get_FA_list(PAL_C16_OCT_C18)),
        'FOH_15 double bond C1-C2': self.SLMs_from_CHEBIs(get_FA_list(FOH_15)),
        'FAS_15': self.SLMs_from_CHEBIs(get_FA_list(FAS_15)),
        'PALOH_C16': self.SLMs_from_CHEBIs(get_FA_list(PALOH_C16)),
        'FAS_79': self.SLMs_from_CHEBIs(get_FA_list(FAS_79)),
        'FAS_85':self.SLMs_from_CHEBIs(get_FA_list(FAS_85))
        }

        # Load FA mapping table
        with importlib.resources.files("swisslipidsreact.package_data").joinpath("FA per class per position.tsv").open("r") as f:
            df = pd.read_csv(f, sep="\t")

        # Adjust for only palmitic acid runs for the test palmitic runs 
        if not curated_fa_list_run:
            for position in positions:
                df.loc[df[position] == 'PAL_C16_OCT_C18', position] = 'PAL_C16'
                df.loc[df[position] == 'FAS_85', position] = 'PAL_C16'
                df.loc[df[position] == 'FAS_15', position] = 'PAL_C16'
                df.loc[df[position] == 'FOH_15', position] = 'PALOH_C16'
                df.loc[df[position] == 'FOH_15 double bond C1-C2', position] = 'PALOH_C16'
                df.loc[df[position] == 'FAS_79', position] = 'PAL_C16'
                df.loc[df[position] == 'SPHINGO_23', position] = 'PAL_C16'

        sn_columns = positions

        # Clean class CHEBI IDs
        df['class_parent_CHEBI'] = df['class_parent_CHEBI'].str.replace("CHEBI:", "").astype(int)

        # Build dictionary only for non-NA components
        chebi_class_to_sn_positions = {
            row['class_parent_CHEBI']: {
                col: row[col] for col in sn_columns
                if pd.notna(row[col]) and row[col] != 'NA'
            }
            for _, row in df.iterrows()
        }

        df_slices = []

        res_filtering = []
        # Main filtering loop
        for chebiid, sn_to_descr in chebi_class_to_sn_positions.items():
            # Here I introduced the flag for degugging to make sure that a certain chebi is 
            # appearing in the dataframe with the descendants
            flag_print = False
            if chebiid == 'test_chebi_id': # Example: 17002
                flag_print = True
            res_chebi_id = [chebiid]
            SLM_classes = self.SLMs_from_CHEBIs([chebiid])
            if flag_print == True:
                print('SLM_classes', SLM_classes)
            if not SLM_classes:
                print(f'No SLM ID found for {chebiid}')
                continue

            df_desc = self.get_lipid_to_descendant_df(SLM_classes)
            if flag_print == True:
                print('descendants before filtering', len(df_desc))
                df_desc.to_csv('test.tsv', sep='\t',index=False)
            res_chebi_id.append(len(df_desc))
            for sn, descr in sn_to_descr.items():
                fa_list = pos_descr_to_FA_list.get(descr)
                if flag_print == True:
                    print('fa list', sn, descr, len(fa_list), 'fatty acids in the list')
                if fa_list is None:
                    print(f'fa list not found for {sn} {descr}')
                    continue  # skip if mapping not found
                df_desc = df_desc[df_desc[sn].isin(fa_list)]
                if flag_print==True:
                    print('descendants after filtering', len(df_desc))
            res_chebi_id.append(len(df_desc))
            res_filtering.append(res_chebi_id)

            if flag_print==True:
                exit()

            # Filter the final_df to retain only rows that are descendants of SLM:000399813 if they have sphing-4-enine base (based on the name)
            glycosphingolipid_slm = "SLM:000399813"

            if glycosphingolipid_slm in SLM_classes:

                # Filter for rows where the rhea_lipid_id corresponds to glycosphingolipid descendants
                descendants_df = df_desc[df_desc['rhea_lipid_id'] == glycosphingolipid_slm]
                # Further filter where Name contains (d18:1(4E)/*)
                df_desc = descendants_df[descendants_df['Name'].str.contains(r'\(d18:1\(4E\)/.+\)', regex=True)]

            df_slices.append(df_desc)
        
        with open(os.path.join(self.output_dir, f'{self.timestamp}_FA_per_class_filtering.tsv'), 'w') as w:
            w.write('Chebiid\tdescendants_total\tdescendants_with_def_components\n')
            for res in res_filtering:
                w.write(f'{res[0]}\t{res[1]}\t{res[2]}\n')

        # Combine all results
        final_df = pd.concat(df_slices, ignore_index=True) if df_slices else pd.DataFrame()

        return final_df
    
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
        df_chebi = self.swisslipids[['Lipid ID', 'CHEBI']].copy()
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