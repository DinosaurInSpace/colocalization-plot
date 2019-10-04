from rdkit import Chem
from rdkit.Chem import PandasTools
from __main__ import * # To access namespace of calling module, e.g. Ipython
import pickle
from pandas import DataFrame as df

# 1. Other half:
# 2. Issue with checking for substrucutre
# 3. dd table of all functional groups...
# 4. Correlate mass spectra with functional groups...

'''
1. Parse log to identfy bad structures
2. Change substruct_variable function and substruct variable to call from config file
3. Source of functional groups:
    https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
    https://www.atmos-chem-phys.net/16/4401/2016/acp-16-4401-2016.pdf
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0225-z
        # Implement algorithm...

'''
# Call as "python molecule_rdkit_pandas_import.py &> log.txt" to write errors to log.

def inchi_smiles(input):

    try:
        molecule = Chem.MolFromInchi(input)
        print("inchi", input '\n')
        output = Chem.MolToSmiles(molecule)
        print("smiles", input, '\n')
        return output

    except:
        print("can't convert", input)
        return 'BrBr'  # Bromine as placeholder for unparcable structures


hmdb_df_static = df(pickle.load(open('hmdb_mols.pickle', 'rb')))

# Make sure original file isn't altered
hmdb_df = hmdb_df_static.copy(deep=True)

# Need to convert to Inchi first!
hmdb_df['Smiles'] = hmdb_df.apply(lambda x: inchi_smiles(x['inchi']), axis=1)

# http://rdkit.org/docs/source/rdkit.Chem.PandasTools.html
# PandasTools.AddMoleculeColumnToFrame(hmdb_df,'Smiles','Molecule')

filename = "hmdb_mol_df.pickle"
outfile = open(filename, "wb")
pickle.dump(hmdb_df, outfile)