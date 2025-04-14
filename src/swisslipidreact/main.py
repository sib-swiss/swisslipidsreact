import pandas as pd
from pyrheadb.RheaDB import RheaDB

from SwissLipids import SwissLipids
from RheaToSwisslipidsDf import RheaToSwisslipidsDf

# list from Fatty tails used the most in glycerolipids.xlsx from Lucila
FAs_str_chebis = """
CHEBI:7896
CHEBI:25629
CHEBI:32372
CHEBI:30823
CHEBI:30245
CHEBI:77222
CHEBI:71589
CHEBI:32395
CHEBI:58562
CHEBI:77016
"""

# palmitate and octadecanoate
FAs_str_chebis = """
CHEBI:7896
CHEBI:25629
"""

# palmitate only
FAs_str_chebis = """
CHEBI:7896
"""

chebis_FAs = [int(i.split(':')[1]) for i in FAs_str_chebis.split() if i]

def split_chebi(chebi):
    """
    Splits a ChEBI identifier (e.g., "CHEBI:1234") into prefix and numeric ID.
    Returns (prefix, numeric ID) or (None, None) on failure.
    """
    try:
        prefix, num = chebi.split(':')
        return pd.Series([prefix, int(num)])
    except Exception:
        return pd.Series([None, None])

def main():
    # ---------- Load and Process RheaDB Data ----------
    rdb = RheaDB()
    rheadf = rdb.rhea_reaction_long_format_smiles_chebi
    rheadf = rheadf[:10000].copy()

    # Remove extract numeric ChEBI IDs
    rheadf[['class', 'chebi_id']] = rheadf['chebiid'].apply(split_chebi)
    rheadf.dropna(subset=['chebi_id'], inplace=True)
    df_rhea_chebi = rheadf[rheadf['class'] != 'POLYMER']
    df_rhea_chebi.loc[:, 'chebi_id'] = df_rhea_chebi['chebi_id'].astype(int)

    sl = SwissLipids()
    sl.read_swisslipids_from_file()

    # Map SwissLipids to ChEBI, and merge with Rhea data
    df_swiss_lipids_chebi = sl.get_only_chebi_sl_df()
    print('# swiss lipids * chebi', len(df_swiss_lipids_chebi))

    FAs = sl.SLMs_from_CHEBIs(chebis_FAs)

    r2sl = RheaToSwisslipidsDf()
    r2sl.get_df_swiss_lipids_chebi_rhea(df_swiss_lipids_chebi, df_rhea_chebi)

    SLMs_in_rhea = r2sl.df_swiss_lipids_chebi_rhea['Lipid ID']

    # Analyse the directed graph of SwissLipid ontology and get all isomeric subspecie per SLM in Rhea
    rhea_lipid_to_descendant_df = sl.get_rhea_lipid_to_descendant_df(SLMs_in_rhea)

    r2sl.get_df_rhea_descendant(rhea_lipid_to_descendant_df, sl.swisslipids)

    r2sl.filter_only_specific_FA_set(rheadf, FA_SLM_list=FAs)

    chebistr = '_'.join([str(i)for i in chebis_FAs])
    r2sl.save_results_per_fa(f'results_{chebistr}.tsv')

main()