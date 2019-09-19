from rdkit import Chem
form rdkit.Chem import PandasTools
from __main__ import *
from contextlib import redirect_stdout
import io
import pickle
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import ast
from rdkit import RDConfig


def structure_converter(raw_structure):
    # Checks for type of encoded structure

    structure = Chem.MolFromInchi(raw_structure)
    if structure == None:
        return Chem.MolFromSmiles('C')
    else:
        return(structure)


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

# Will this not work???  Type of series is wrong...
hmdb_df = df(pickle.load(open('hmdb_mols.pickle', 'rb')))

print(type(hmdb_df['inchi']))

hmdb_df['rd_struc'] = hmdb_df.apply(structure_converter(hmdb_df['inchi']))
'''
substructs = ["-OH", "-NH2", "-CO2H"]

for target_name in substructs:
    hmdb_df[target_name] = substruct_target(target_name)

    present = target_name + "_present"
    hmdb_df[present] = hmdb_df.apply(lambda x: nl_finder(x['rd_struct'], x[target_name]), axis=1)
'''
filename = "hmdb_df.pickle"
outfile = open(filename, "wb")
pickle.dump(hmdb_df, outfile)
