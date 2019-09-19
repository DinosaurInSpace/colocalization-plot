from rdkit import Chem
from rdkit.Chem import PandasTools
from __main__ import * # To access namespace of calling module, e.g. Ipython
from contextlib import redirect_stdout
import io
import pickle
from pandas import DataFrame as df


# 1. Save log of when structures fail, and triage
# 2. Change substruct_variable function and substruct variable to call from config file
# 3. dd table of all functional groups...
# 4. Correlate mass spectra with functional groups...

def substruct_target(target_name):
    # Smarts is more specific for searching

    # Hydroxyl in alcohol
    hydroxy = Chem.MolFromSmarts('[CX4][OX2H]')

    # Primary or secondary amine, not amide.
    amine = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')  # Similiar issue, especially oxidation.

    # Carboxylic acid or conjugate base.
    carboxyl = Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]')

    targets = {"-OH": hydroxy, "-NH2": amine, "-CO2H": carboxyl}

    return targets[target_name]


def nl_finder(structure, target):
    # Can find arbitrary structures
    # Has SubstructMatch only delivers to standard out, instead capture as file.
    with io.StringIO() as buf, redirect_stdout(buf):
        output = buf.getvalue()
        print(structure.HasSubstructMatch(target))
        return output


def inchi_smiles(input):

    try:
        molecule = Chem.MolFromInchi(input)
        output = Chem.MolToSmiles(molecule)
        return output
    except:
        return 'C'  # Methane as placeholder for unparcable structures


hmdb_df_static = df(pickle.load(open('hmdb_mols.pickle', 'rb')))
hmdb_df = hmdb_df_static.copy(deep=True)

#Need to convert to Inchi first!
hmdb_df['Smiles'] = hmdb_df.apply(lambda x: inchi_smiles(x['inchi']), axis=1)

PandasTools.AddMoleculeColumnToFrame(hmdb_df,'Smiles','Molecule')


substructs = ["-OH", "-NH2", "-CO2H"]

for target_name in substructs:
    hmdb_df[target_name] = substruct_target(target_name)

    present = target_name + "_present"
    hmdb_df[present] = hmdb_df.apply(lambda x: nl_finder(x['Molecule'], x[target_name]), axis=1)

filename = "hmdb_df.pickle"
outfile = open(filename, "wb")
pickle.dump(hmdb_df, outfile)