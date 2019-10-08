from contextlib import redirect_stdout
import io
from pandas import DataFrame as df
import pickle
from rdkit import Chem
from rdkit.Chem import PandasTools
from structures_to_search_n5 import target_structures #n5

# Run command as: "python hmdb_structure_parser.py 2>&1 | tee log.txt" to get log and stdout
# 3.5 minutes run time

def open_pickle(name):
    current_df_static = df(pickle.load(open(name, 'rb')))
    current_df_open = current_df_static.copy(deep=True)
    return current_df_open


def inchi_smiles(input_mol):
    # Converts molecules from INCHI to SMILES

    try:
        print('MolFromInchi', input_mol)
        molecule = Chem.MolFromInchi(input_mol)
        print('MolToSmiles', input_mol)
        output_smiles = Chem.MolToSmiles(molecule)
        return output_smiles

    except:
        print('ERROR')
        return '[Si]'  # Silicon as placeholder for unparcable structures


def substruct_target(active_df):
    # Adds target structures as rd_object to dataframe

    for target, substruct in target_structures.items():
            substruct_object = Chem.MolFromSmarts(substruct)
            target_name = str(target) + '_target'
            active_df[target_name] = substruct_object
    return active_df


def pandas_structure(active_df):
    # Converts INCHI input file to smiles, then adds rd_object to dataframe'''

    active_df['Smiles'] = active_df.apply(lambda x: inchi_smiles(x['inchi']), axis=1)
    PandasTools.AddMoleculeColumnToFrame(active_df, 'Smiles', 'Molecule')
    return active_df


def pickle_out(out_df, outname):
    outfile = open(outname, "wb")
    pickle.dump(out_df, outfile)
    outfile.close()


current_df = open_pickle('hmdb_mols.pickle')
current_df = substruct_target(current_df)
current_df = pandas_structure(current_df)

filename = "hmdb_df_n5.pickle" # Out name editied
pickle_out(current_df, filename)