#%% Setup
import pickle
import pandas as pd
TAXONOMY_FIELDS = ['kingdom', 'super_class', 'class', 'sub_class', 'direct_parent', 'molecular_framework']

#%% Parse HMDB & dump to hmdb_mols.pickle
from xml.etree.cElementTree import iterparse

def get_one_text(elem, q):
    node = elem.find(q)
    return node.text if node is not None else None

def get_many_text(elem, q):
    return [node.text for node in elem.findall(q)]

NS = '{http://www.hmdb.ca}'

root = None
mols = []
for event, elem in iterparse('hmdb_metabolites.xml', events=('start','end')):
    if root is None:
        root = elem
    if event == 'end' and elem.tag == f'{NS}metabolite':
        fields = {
            'id': get_one_text(elem, f'{NS}accession'),
            'formula': get_one_text(elem, f'{NS}chemical_formula'),
            'inchi': get_one_text(elem, f'{NS}inchi'),
            'mol_name': get_one_text(elem, f'{NS}name'),
            **dict((field, get_one_text(elem, f'{NS}taxonomy/{NS}' + field)) for field in TAXONOMY_FIELDS),
            'substituents': get_many_text(elem, f'{NS}taxonomy/{NS}substituents/{NS}substituent'),
            'cellular_locations': get_many_text(elem, f'{NS}biological_properties/{NS}cellular_locations/{NS}cellular'),
            'biospecimen_locations': get_many_text(elem, f'{NS}biological_properties/{NS}biospecimen_locations/{NS}biospecimen'),
            'tissue_locations': get_many_text(elem, f'{NS}biological_properties/{NS}tissue_locations/{NS}tissue'),
            'pathways': get_many_text(elem, f'{NS}biological_properties/{NS}pathways/{NS}pathway/{NS}name'),
        }
        mols.append(fields)
        if len(mols) % 1000 == 0: print(len(mols))

        elem.clear()
        root.clear()

# Dump the output to disk so it can be loaded more quickly later
pickle.dump(mols, open('hmdb_mols.pickle', 'wb'))

#%%
# Reload
mols = pickle.load(open('hmdb_mols.pickle', 'rb'))