from pandas import DataFrame as df
import pickle
from rdkit import Chem
from rdkit.Chem import Fragments

# from rdkit.Chem import PandasTools

"""
Status:
1. Executes 2 min. , runs without error...
2. To do: delete columns ending in "_target"
3. Need to use sanitized file re: line zero
4. Add input file and output file names at CL
4. Add search objects to DF here first instead of hmdb structure parser
5. Add search for none type objects in "structures to search"
5. To do: Correlate mass spectra with functional groups...
"""


def header_targets(y_df):
    columns = y_df.columns
    header_list = []

    for header in columns:
        if '_target' in header:
            header_list.append(header)
        else:
            continue
    return header_list


def nl_finder(structure, target):
    # Can find arbitrary structures

    search_result = structure.HasSubstructMatch(target)
    return search_result


def struct_pandas_search(x_df, headers):
    for head in headers:
        res = head + "_present"
        x_df[res] = x_df.apply(lambda x: nl_finder(x['Molecule'], x[head]), axis=1)
    return x_df


current_df = df(pickle.load(open('hmdb_df_n5.pickle', 'rb'))) # Name changed
current_headers = header_targets(current_df)
final_df = struct_pandas_search(current_df, current_headers)

# Can only find specific hardcoded structures in RDkit
# final_df['CA_alaph_alt'] = final_df.apply(lambda x: bool(Chem.Fragments.fr_Al_COO(x['Molecule'])), axis=1)
# final_df['OH_alaph_alt'] = final_df.apply(lambda x: bool(Chem.Fragments.fr_Al_OH(x['Molecule'])), axis=1)

# OCD:
final_df = final_df.reindex(sorted(final_df.columns), axis=1)

outname = "hmdb_df_n5_searched.pickle" # Name changed
outfile = open(outname, "wb")
pickle.dump(final_df, outfile)
outfile.close()
