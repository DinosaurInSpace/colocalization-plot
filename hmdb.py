#%% Setup
import pickle
import pandas as pd
TAXONOMY_FIELDS = ['kingdom', 'super_class', 'class', 'sub_class', 'direct_parent', 'molecular_framework']
LIPID_CLASS = 'Lipids and lipid-like molecules'

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
for event, elem in iterparse('/home/lachlan/Downloads/HMDB/hmdb_metabolites.xml', events=('start','end')):
    if root is None:
        root = elem
    if event == 'end' and elem.tag == f'{NS}metabolite':
        fields = {
            'id': get_one_text(elem, f'{NS}accession'),
            'formula': get_one_text(elem, f'{NS}chemical_formula'),
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

pickle.dump(mols, open('/home/lachlan/dev/notebooks/metaspace-mol-cloud/hmdb_mols.pickle', 'wb'))

#%%
mols = pickle.load(open('/home/lachlan/dev/notebooks/metaspace-mol-cloud/hmdb_mols.pickle', 'rb'))

#%%
basic_df = pd.DataFrame(mols)[['id', 'formula', 'mol_name', *TAXONOMY_FIELDS]].set_index('id')
cloc_df = pd.DataFrame([{'id': m['id'], **dict((l, 1) for l in m['cellular_locations'])} for m in mols]).set_index('id')
bloc_df = pd.DataFrame([{'id': m['id'], **dict((l, 1) for l in m['biospecimen_locations'])} for m in mols]).set_index('id')
tloc_df = pd.DataFrame([{'id': m['id'], **dict((l, 1) for l in m['tissue_locations'])} for m in mols]).set_index('id')
#%%
basic_stats = basic_df.describe()
taxonomy_freqs = dict((k, basic_df[k].value_counts().sort_values(ascending=False)) for k in TAXONOMY_FIELDS)
cloc_stats = cloc_df.describe()
bloc_stats = bloc_df.describe()
tloc_stats = tloc_df.describe()

#%%
bloc_by_formula = basic_df[['formula']].join(bloc_df).set_index('formula').groupby(level=0).sum()
bloc_by_formula = bloc_by_formula.applymap(lambda val: val and 1 or 0)
tloc_by_formula = basic_df[['formula']].join(bloc_df).set_index('formula').groupby(level=0).sum()
tloc_by_formula = tloc_by_formula.applymap(lambda val: val and 1 or 0)
#%%
def tag_lipids(df):
    if (df.super_class == LIPID_CLASS).any():
        is_lipid = 'Lipid' if (df.super_class == LIPID_CLASS).all() else 'Maybe lipid'
        lipid_type = df[df.super_class == LIPID_CLASS]['class'].value_counts().index[0]
        if lipid_type not in ['Glycerolipids','Glycerophospholipids','Fatty Acyls','Prenol lipids','Steroids and steroid derivatives']:
            lipid_type = 'Other lipids'
    else:
        is_lipid = 'Non-lipid'
        lipid_type = 'Non-lipid'

    def most_common(col):
        c = df[col].value_counts()
        if len(c) > 0:
            return c.index[0]
    return pd.DataFrame({'formula': [df.formula.iloc[0]],
                         'is_lipid': [is_lipid],
                         'lipid_type': [lipid_type],
                         'super_class': most_common('super_class'),
                         'class': most_common('class'),
                         'sub_class': most_common('sub_class')})
lipids_df = basic_df.groupby('formula', group_keys=False).apply(tag_lipids)
#%%

ion_formulas = pd.concat([pd.DataFrame([[formula, (formula + adduct)] for adduct in ['+H','+Na','+K','-H','+Cl']], columns=['formula','ion']) for formula in basic_df.formula.drop_duplicates()])
#%%
basic_by_ion = ion_formulas.merge(basic_df, on='formula').set_index('ion')
lipids_by_ion = ion_formulas.merge(lipids_df, on='formula').set_index('ion')
# bloc_by_ion = bloc_by_formula.merge(ion_formulas, left_index=True, right_on='formula').set_index('ion').drop(columns=['formula'])
# tloc_by_ion = tloc_by_formula.merge(ion_formulas, left_index=True, right_on='formula').set_index('ion').drop(columns=['formula'])
#%%
basic_by_ion.to_csv('/home/lachlan/dev/notebooks/metaspace-mol-cloud/hmdb_by_ion.csv')
lipids_df.to_csv('/home/lachlan/dev/notebooks/metaspace-mol-cloud/lipids_by_formula.csv', index=False)
lipids_by_ion.to_csv('/home/lachlan/dev/notebooks/metaspace-mol-cloud/lipids_by_ion.csv')
bloc_by_ion.to_csv('/home/lachlan/dev/notebooks/metaspace-mol-cloud/hmdb_biospecimen_loc.csv')
tloc_by_ion.to_csv('/home/lachlan/dev/notebooks/metaspace-mol-cloud/hmdb_tissue_loc.csv')
#%% Make hmdb_with_array_fields
def merge_list_fields(group):
    return pd.DataFrame([{
        'formula': group.name,
        **dict((k, '|'.join(sorted(set(group[k]) - {None}))) for k in TAXONOMY_FIELDS),
        **dict((k, '|'.join(sorted(set(i for sublist in group[k] for i in sublist if i)))) for k in ['cellular_locations', 'biospecimen_locations', 'tissue_locations']),
    }], index='formula')
flat_array_df = (pd.DataFrame(mols)
                 .groupby('formula')
                 .apply(merge_list_fields)
                 .merge(ion_formulas, left_index=True, right_on='formula')
                 .set_index('ion', drop=True))
flat_array_df.to_csv('/home/lachlan/dev/notebooks/metaspace-mol-cloud/hmdb_with_array_fields.csv')
#%%

#%%

#%%

#%%
