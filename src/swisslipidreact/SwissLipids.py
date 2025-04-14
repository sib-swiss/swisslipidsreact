import pandas as pd
import os
import networkx as nx

SLMs_in_rhea = []

class SwissLipids():

    def __init__(self):
        pass

    # ---------- Load SwissLipids Data ----------
    def read_swisslipids_from_file(self):
        self.swisslipids = pd.read_csv(
            os.path.join('..','..', 'data','lipids.tsv'), sep='\t', encoding='latin-1',
            usecols=['Lipid ID', 'CHEBI', 'Level', 'Lipid class*', 'Components*', 'SMILES (pH7.3)'],
            dtype={'Lipid ID': str, 'CHEBI': str, 'Level': str, 'Lipid class*': str,
                'Components*': str, 'SMILES (pH7.3)': str}
        )

    def get_only_chebi_sl_df(self):
        """
        Extracts and processes CHEBI identifiers from the SwissLipids dataframe.
        Returns a dataframe mapping 'Lipid ID' to individual numeric ChEBI IDs.
        """
        df_chebi = self.swisslipids[['Lipid ID', 'CHEBI']].copy()
        df_chebi.dropna(subset=['CHEBI'], inplace=True)
        df_chebi['CHEBI'] = df_chebi['CHEBI'].str.split('|')
        df_chebi = df_chebi.explode('CHEBI')
        df_chebi.dropna(subset=['CHEBI'], inplace=True)
        df_chebi['chebi_id'] = df_chebi['CHEBI'].apply(lambda x: int(x.replace('CHEBI:', '')))
        self.df_chebi = df_chebi
        return df_chebi
   
    def SLMs_from_CHEBIs(self, list_of_chebi_ids):
        return self.df_chebi[self.df_chebi['chebi_id'].isin(list_of_chebi_ids)]['Lipid ID'].to_list()
    
    def get_isomeric_subspecies_table(self):
        # Filter isomeric subspecies
        isomeric_subspecies = self.swisslipids[self.swisslipids['Level'] == 'Isomeric subspecies']['Lipid ID'].tolist()
        self.df_isomeric_subspecies = pd.DataFrame(isomeric_subspecies, columns=['Lipid ID'])

    # ---------- Lipid Class Graph Analysis ----------
    def get_rhea_lipid_to_descendant_df(self, SLMs_in_rhea):
        # Build a directed graph of lipid class relationships
        G_lipid_class = nx.from_pandas_edgelist(
            self.swisslipids, source='Lipid class*', target='Lipid ID', create_using=nx.DiGraph()
        )

        # Find descendants for each lipid in the Rhea-SwissLipids merged set
        list_id_descendant = []
        for lipid_id in set(SLMs_in_rhea):
            descendants = nx.descendants(G_lipid_class, lipid_id)
            list_id_descendant.extend((lipid_id, i) for i in descendants)

        # Create dataframe of isomeric subspecies relationships
        rhea_lipid_to_descendant_df = pd.DataFrame(list_id_descendant, columns=[
            'rhea_lipid_id', 'isomeric_subspecies_descendant_lipid_id'
        ])
        print('all descendants of swiss lipids * Rhea', len(rhea_lipid_to_descendant_df))
        print('all descendants of swiss lipids in Rhea', len(set(rhea_lipid_to_descendant_df['isomeric_subspecies_descendant_lipid_id'])))

        # Filter to retain only isomeric subspecies
        self.get_isomeric_subspecies_table()

        rhea_lipid_to_descendant_df = rhea_lipid_to_descendant_df.merge(
            self.df_isomeric_subspecies, left_on='isomeric_subspecies_descendant_lipid_id',
            right_on='Lipid ID', how='inner'
        )

        rhea_lipid_to_descendant_df.drop(columns=['Lipid ID'], inplace=True)
        print('isomeric subspecies descendants of lipids * Rhea', len(rhea_lipid_to_descendant_df))
        print('isomeric subspecies descendants of lipids in Rhea', len(set(rhea_lipid_to_descendant_df['isomeric_subspecies_descendant_lipid_id'])))

        return rhea_lipid_to_descendant_df
