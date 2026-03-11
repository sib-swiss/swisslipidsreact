import logging
import pandas as pd
import os
import networkx as nx
import importlib.resources
from platformdirs import user_cache_dir

from .utils import DEBUG, debug_df_first_row
from .FA_lists import positions, get_FA_list, FAS_15, FAS_85, FAS_79, FOH_15, PAL_C16, PALOH_C16, PAL_C16_OCT_C18, SPHINGO_23

logger = logging.getLogger(__name__)

flag_fast_exec = True

class SwissLipids():

    def __init__(self, output_dir='', timestamp='now'):
        # CACHE: Create a project-specific subdirectory in the user's cache directory.
        self.cache_dir = user_cache_dir("swisslipidsreact")
        os.makedirs(self.cache_dir, exist_ok=True)
        self.output_dir = output_dir
        self.timestamp=timestamp

    def pre_process_swisslipids(self, output_file):
        """
        This function adds (free fatty acid) as component to FA and FA-COA and n-acyls to acylethanolamines
        """
        logger.info( "Pre-processing SwissLipids file" )
        
        with importlib.resources.files("swisslipidsreact.package_data").joinpath("lipids.tsv").open("rb") as f:
            self.swisslipids = pd.read_csv(f, sep="\t", encoding="latin-1",
                usecols=['Name', 'Lipid ID', 'CHEBI', 'Level', 'Lipid class*', 'Components*', 'SMILES (pH7.3)'],
                dtype={'Lipid ID': str, 'CHEBI': str, 'Level': str, 'Lipid class*': str,
                       'Components*': str, 'SMILES (pH7.3)': str})

        # Assign level 'Isomeric subspecies' for all lipids that are not class lipids (compounds with *).
        self.swisslipids.loc[~self.swisslipids['SMILES (pH7.3)'].isna() & 
        self.swisslipids['Level'].isna() & ~self.swisslipids['SMILES (pH7.3)'].str.contains('*', regex=False, na=False), 'Level'] = 'Isomeric subspecies'

        # Normalize prime characters to ASCII apostrophe (')
        prime_variants = ['′', 'ʹ', '´']  # U+2032, U+02B9, U+00B4
        for variant in prime_variants:
            self.swisslipids['Components*'] = self.swisslipids['Components*'].str.replace(variant, "'", regex=False)

        # Helper function to clean the Components* column
        def __clean_components(row):
            """
            Removes class from the components
            """
            components = row['Components*']
            lipid_class = row['Lipid class*']
            
            if pd.isna(components):
                return components  # Leave NaN as is

            # Split and filter.
            ### TODO: "BUG" - 119 rows have multiple lipid_class, e.g. "SLM:000000261 | SLM:000056434"
            ### -> There seem to be no cases in lipids.tsv, why is she doing this?
            ### Lipid class* that matches Components* (estimation, not splitting the 119 cases of Lipid class* with '|' )
            ### Ex: Isomeric subspecies	SLM:000000347	SLM:000001101 (sn1) / SLM:000001200 (sn2)
            ### tail -n+2 lipids.tsv | perl -ne 'chomp; @a=split(/\t/); print "$a[1]\t$a[5]\t$a[7]\n" if $a[7] and ($a[5] =~ m|$a[7]|)' | wc -l
            ### 0
            filtered = [comp.strip() for comp in components.split(' / ') if lipid_class not in comp]
            return ' / '.join(filtered) if filtered else None

        # Apply __clean_components function
        self.swisslipids['Components*'] = self.swisslipids.apply(__clean_components, axis=1)

        self.get_isomeric_subspecies_table() # generates self.df_isomeric_subspecies
        
        self.get_lipid_class_graph()

        # Assign free fatty acid as component for all FA and FA-CoA
        FA = "SLM:000000984"
        F_alcohol = "SLM:000390053"

        FAs = self.build_df_slm_class2iso([FA])['Lipid ID'].tolist()
        F_alcohols = self.build_df_slm_class2iso([F_alcohol])['Lipid ID'].tolist()
        
        # Assign (free fatty acid) to Components* for FAs
        self.swisslipids.loc[self.swisslipids['Lipid ID'].isin(FAs), 'Components*'] = (
        self.swisslipids['Lipid ID'] + ' (free fatty acid)' )

        # Assign (free fatty alcohol) to Components* for FAs
        self.swisslipids.loc[self.swisslipids['Lipid ID'].isin(F_alcohols), 'Components*'] = (
        self.swisslipids['Lipid ID'] + ' (free fatty alcohol)' )

        # Helper function for readability.
        def __load_package_tsv(module, filename, sep="\t"):
            with importlib.resources.files(module).joinpath(filename).open("r") as f:
                return pd.read_csv(f, sep=sep)

        # Now load the files:
        df_facoa_to_ffa = __load_package_tsv(
            "swisslipidsreact.package_data.components_correction", 
            "free fatty acids per compound.tsv"
        )
        comp_to_ffa_dict = dict(zip(df_facoa_to_ffa['Lipid ID'], df_facoa_to_ffa['free fatty acid']))

        df_comp_to_nacyl_dict = __load_package_tsv(
            "swisslipidsreact.package_data.components_correction",
            "nacyls per compound.tsv"
        )
        comp_to_nacyl_dict = dict(zip(df_comp_to_nacyl_dict['Lipid ID'], df_comp_to_nacyl_dict['n-acyl']))

        df_comp_to_sn1_dict = __load_package_tsv(
            "swisslipidsreact.package_data.components_correction",
            "sn1 per compound.tsv"
        )
        comp_to_sn1_dict = dict(zip(df_comp_to_sn1_dict['Lipid ID'], df_comp_to_sn1_dict['sn1']))

        df_comp_to_sn2_dict = __load_package_tsv(
            "swisslipidsreact.package_data.components_correction",
            "sn2 per compound.tsv"
        )
        comp_to_sn2_dict = dict(zip(df_comp_to_sn2_dict['Lipid ID'], df_comp_to_sn2_dict['sn2']))
        
        def __assign_components(row, comp_to_ffa_dict, comp_to_nacyl_dict, comp_to_sn1_dict, comp_to_sn2_dict):
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

        # Assign values to Components* for CoAs.
        # NB: This code adds components to lipids that were not enumerated by us, but that exist in ChEBI.
        # TODO: This may not be 100% correct, we need to find a Rhea for it.
        #       Alan thought the code just ignored the things we have not enumerated.
        self.swisslipids['Components*'] = self.swisslipids.apply(__assign_components, axis=1, args = [comp_to_ffa_dict, comp_to_nacyl_dict, comp_to_sn1_dict, comp_to_sn2_dict, ])

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

        self.swisslipids.to_csv(output_file, sep="\t", index=False)
    
    # ---------- Load SwissLipids Data ----------
    def read_swisslipids_from_file(self):

        # TODO: This is fast, does not need caching?
        # TODO: Remove hard-coded flag flag_fast_exec
        f_cache_lipids_preprocessed = os.path.join(self.cache_dir, 'lipids_preprocessed.tsv')
        if not os.path.exists(os.path.join(f_cache_lipids_preprocessed)):
            logger.info( "Generating pre-processed SwissLipids file" )
            self.pre_process_swisslipids(f_cache_lipids_preprocessed)

        # TODO: This is fast, does not need caching?
        lipids_components_split_cache_file = os.path.join(self.cache_dir, 'lipids_components_split.tsv')
        if not os.path.exists(lipids_components_split_cache_file) or flag_fast_exec==False:
            logger.info( "Extracting components into positions" )
            self.swisslipids = pd.read_csv(
                os.path.join(self.cache_dir, "lipids_preprocessed.tsv"), sep='\t', encoding='latin-1',
                usecols=['Name', 'Lipid ID', 'CHEBI', 'Level', 'Lipid class*', 'Components*', 'SMILES (pH7.3)'],
                dtype={'Lipid ID': str, 'CHEBI': str, 'Level': str, 'Lipid class*': str,
                       'Components*': str, 'SMILES (pH7.3)': str}
            )

            # TODO: This is duplicated in pre_process_swisslipids?
            self.swisslipids.loc[~self.swisslipids['SMILES (pH7.3)'].isna() & 
            self.swisslipids['Level'].isna() & ~self.swisslipids['SMILES (pH7.3)'].str.contains('*', regex=False, na=False), 'Level'] = 'Isomeric subspecies'

            # TODO: Ask Lucila what this is
            # malonyl-CoA(5−)	CHEBI:57384
            # acetyl-CoA(4−)	CHEBI:57288
            # glycerone phosphate(2−)	CHEBI:57642
            # (S)-methylmalonyl-CoA(5−)	CHEBI:57327
            # acetoacetate	CHEBI:13705
            # dicarboxylic acid dianion CHEBI:28965
            self.swisslipids=self.swisslipids[~self.swisslipids['CHEBI'].isin(['57384', '57288', '57642', '57327', '13705', '28965'])]

            # Define possible positions
            # Initialize columns to None
            for pos in positions:
                self.swisslipids[pos] = None

            # Work only on rows where Level == 'Isomeric subspecies'
            is_iso = self.swisslipids['Level'] == 'Isomeric subspecies'
            df_iso = self.swisslipids[is_iso].copy()

            # The 'Components*' column contains 1-n components, each with 1-m positions :-(
            # Example: SLM:000001124 (sn2) / SLM:000001208 (sn1 or sn1')
            # Here we create one column for each position and store the corresponding component there.
            pattern = r'(SLM:\d+)\s+\(([^)]+)\)'
            extracted = df_iso['Components*'].str.extractall(pattern)
            extracted.columns = ['SLM', 'positions']
            extracted['row'] = extracted.index.get_level_values(0)

            # Expand positions string (e.g. sn1 or sn1').
            extracted = extracted.assign(position=extracted['positions'].str.split(r'\s*(?:or|,|and)\s*')).explode('position')
            extracted['position'] = extracted['position'].str.strip() # Correct whitespace errors.

            # Keep only the defined positions (there should be no others).
            extracted = extracted[extracted['position'].isin(positions)]
            if len(extracted[~extracted['position'].isin(positions)] > 0):
                logger.warning("Unknown positions: %s\n", extracted[~extracted['position'].isin(positions)])

            # Drop duplicates to enforce a single SLM per position.
            extracted = extracted.drop_duplicates(subset=['row', 'position'])

            # Pivot to wide format
            components_wide = extracted.pivot(index='row', columns='position', values='SLM')

            # Add the new columns to the original df (only for Isomeric subspecies rows).
            for pos in positions:
                self.swisslipids.loc[is_iso, pos] = df_iso.index.map(components_wide.get(pos))
            
            self.swisslipids.to_csv(lipids_components_split_cache_file, sep='\t', index=False)
            
        else:
            logger.info( "Reading pre-processed SwissLipids file with components in positions" )
            self.swisslipids = pd.read_csv(lipids_components_split_cache_file, sep='\t', low_memory=False)
        
        # Filter to retain only isomeric subspecies in a seeparate df
        # TODO: This is duplicated in pre_process_swisslipids?
        self.get_isomeric_subspecies_table() # generates self.df_isomeric_subspecies
        
        self.get_lipid_class_graph()
    
    def get_lipid_class_graph(self):
        """
        Builds a directed graph of lipid class relationships.
        """
        # 'Lipid class*' may contain a pipe-separated list of values: Expand each value into a separate row.
        df_expanded = self.swisslipids.assign(
            **{'Lipid class*': self.swisslipids['Lipid class*'].str.split('|')}
        ).explode('Lipid class*')
        # Strip extra whitespace (bug in lipids.tsv).
        df_expanded['Lipid class*'] = df_expanded['Lipid class*'].str.strip()

        # Now build the graph.
        self.G_lipid_class = nx.from_pandas_edgelist(
            df_expanded, source='Lipid class*', target='Lipid ID', create_using=nx.DiGraph()
        )
    

    def filter_curated_biologically_relevant_isomeric_subspecies_only(self, curated_fa_list_run=True):
        """
        Creates df prefiltered to only have lipids specific for a particular positions / FA
        Columns:
        ['class_slm_id', 'iso_slm_id', 'Lipid ID',
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
        all_lipids_considered = []
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

            df_desc = self.build_df_slm_class2iso(SLM_classes)
            all_lipids_considered.extend(df_desc['Lipid ID'].tolist())
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

                # Filter for rows where the class_slm_id corresponds to glycosphingolipid descendants
                descendants_df = df_desc[df_desc['class_slm_id'] == glycosphingolipid_slm]
                # Further filter where Name contains (d18:1(4E)/*)
                df_desc = descendants_df[descendants_df['Name'].str.contains(r'\(d18:1\(4E\)/.+\)', regex=True)]

            df_slices.append(df_desc)
        
        with open(os.path.join(self.output_dir, f'{self.timestamp}_FA_per_class_filtering.tsv'), 'w') as w:
            w.write('Chebiid\tdescendants_total\tdescendants_with_def_components\n')
            for res in res_filtering:
                w.write(f'{res[0]}\t{res[1]}\t{res[2]}\n')

        # Combine all results
        final_df = pd.concat(df_slices, ignore_index=True) if df_slices else pd.DataFrame()

        return final_df, all_lipids_considered
    
    # ---------- Lipid Class Graph Analysis ----------
    def build_df_slm_class2iso(self, class_SLMs):

        # Find the descendants of each provided class SLM.
        list_id_descendant = []
        for slm_id in set(class_SLMs):
            descendants = nx.descendants(self.G_lipid_class, slm_id)
            list_id_descendant.extend((slm_id, i) for i in descendants)

        # Create dataframe of class to isomeric subspecies relationships.
        df_slm_class2iso = pd.DataFrame(list_id_descendant, columns=[
            'class_slm_id', 'iso_slm_id'
        ])

        df_slm_class2iso = df_slm_class2iso.merge(
            self.df_isomeric_subspecies, left_on='iso_slm_id',
            right_on='Lipid ID', how='inner'
        )
        return df_slm_class2iso
    
    def build_df_slm2chebi(self):
        """
        Builds a map of SLM IDs to numeric ChEBI IDs.
        Deals with pipe-separated 'CHEBI' values and erroneous prefixes.
        """
        df_slm2chebi = self.swisslipids[['Lipid ID', 'CHEBI']].copy()
        df_slm2chebi.dropna(subset=['CHEBI'], inplace=True)
        # 'CHEBI' may contain a pipe-separated list of values: Expand each value into a separate row.
        # NB: There are only 3 rows: 2 are likely errors due to CHEBI merges, the 3rd is a "duplication":
        # 74546 | 82922
        # 82731 | CHEBI:82731
        # 17336 | 83228
        df_slm2chebi['CHEBI'] = df_slm2chebi['CHEBI'].astype(str).str.split('|')
        df_slm2chebi = df_slm2chebi.explode('CHEBI')
        df_slm2chebi.dropna(subset=['CHEBI'], inplace=True)
        # Bug in lipids.tsv: There is one row with CHEBI:82731 instead of 82731.
        df_slm2chebi['chebi_id'] = df_slm2chebi['CHEBI'].apply(lambda x: int(float(x.replace('CHEBI:', ''))) if isinstance(x, str) else int(x))
        self.df_slm2chebi = df_slm2chebi
        return df_slm2chebi
   
    def SLMs_from_CHEBIs(self, list_of_chebi_ids):

        return self.df_slm2chebi[self.df_slm2chebi['chebi_id'].isin(list_of_chebi_ids)]['Lipid ID'].tolist()
    
    def get_isomeric_subspecies_table(self):
        # Filter isomeric subspecies
        isomeric_subspecies = self.swisslipids[self.swisslipids['Level'] == 'Isomeric subspecies']
        self.df_isomeric_subspecies = isomeric_subspecies.copy()
        logger.info( "LENGTH: %8d df_isomeric_subspecies", len(self.df_isomeric_subspecies) )
        logger.debug( "FORMAT: df_isomeric_subspecies\n%s", debug_df_first_row(self.df_isomeric_subspecies) )
        if DEBUG > 1:
            self.df_isomeric_subspecies.to_csv(os.path.join(self.output_dir, 'DEBUG_df_isomeric_subspecies.tsv'), sep="\t", header=True, index=False)
