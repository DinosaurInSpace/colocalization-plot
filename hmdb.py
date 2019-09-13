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


#%% Build dataframes of features for each molecule
LIPID_CLASS = 'Lipids and lipid-like molecules'
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
def mols_to_ratios(group):
    nodes = [
        group.super_class,
        group.super_class + '/' + group['class'],
        group.super_class + '/' + group['class'] + '/' + group.sub_class,
    ]
    s = (pd.concat([node.value_counts() / len(group) for node in nodes])
         .apply(lambda val: 'Likely' if val == 1 else 'Possibly' if val != 0 else 'Not')
         .astype('category')
         .rename(group.name))
    return pd.DataFrame({'tree': s.index, 'ratio': s.values})

ratio_array_df = (pd.DataFrame(mols)
                 .groupby('formula')
                 .apply(mols_to_ratios))
ratio_array_df = ratio_array_df.reset_index().pivot(values='ratio', index='formula', columns='tree')
ratio_array_df = ratio_array_df.merge(ion_formulas, left_index=True, right_on='formula')
ratio_array_df = ratio_array_df.set_index('ion', drop=True)
ratio_array_df.to_pickle('/home/lachlan/dev/notebooks/metaspace-mol-cloud/hmdb_possibility_pivot.pickle')
#%%
def all_or_sum(series):
    return 'Likely' if series.all() else 'Possibly' if series.any() else 'Not'
def classes_of_interest(group):
    s = {
        'Amino acids, peptides, and analogues': group.sub_class == 'Amino acids, peptides, and analogues',
        'Fatty acids and conjugates': group.sub_class == 'Fatty acids and conjugates',
        'Long-chain fatty acids': group.direct_parent == 'Long-chain fatty acids',
        'Lysophospholipids': group.direct_parent.isin({
            "1-acyl-sn-glycero-3-phosphoethanolamines",
            "2-acyl-sn-glycero-3-phosphoethanolamines",
            "1-(1Z-alkenyl)-glycero-3-phosphoethanolamines",
            "1-acyl-sn-glycero-3-phosphocholines",
            "1-(1Z-alkenyl)-glycero-3-phosphocholines",
            "1-acylglycerol-3-phosphates",
            "2-acylglycerol-3-phosphates",
            "1-acyl-sn-glycerol-3-phosphoinositols",
        }),
        'Ceramides': group.sub_class == 'Ceramides',
        'Phosphosphingolipids': group.sub_class == 'Phosphosphingolipids',
        'Sphingolipids': group['class'] == 'Sphingolipids',
        'Glycerophospholipids': group['class'] == 'Glycerophospholipids',
        'Glycerophosphocholines': group.sub_class == 'Glycerophosphocholines',
        'Glycerophosphoethanolamines': group.sub_class == 'Glycerophosphoethanolamines',
        'Glycerophosphates': group.sub_class == 'Glycerophosphates',
        'Glycerophosphoglycerols': group.sub_class == 'Glycerophosphoglycerols',
        'Glycerophosphoinositols': group.sub_class == 'Glycerophosphoinositols',
        'Glycerolipids': group['class'] == 'Glycerolipids',
        'Triradylcglycerols': group.sub_class == 'Triradylcglycerols',
        'Diradylglycerols': group.sub_class == 'Diradylglycerols',
        'Monoradylglycerols': group.sub_class == 'Monoradylglycerols',
    }
    s = pd.Series(dict([(k, all_or_sum(v)) for k,v in s.items()]))
    return pd.DataFrame({'tree': s.index, 'ratio': s.values})

class_array_df = (pd.DataFrame(mols)
                 .groupby('formula')
                 .apply(classes_of_interest))
class_array_df = class_array_df.reset_index().pivot(values='ratio', index='formula', columns='tree')
class_array_df = class_array_df.merge(ion_formulas, left_index=True, right_on='formula')
class_array_df = class_array_df.set_index('ion', drop=True)
class_array_df = class_array_df.drop(columns='formula')
class_array_df.to_csv('/home/lachlan/dev/notebooks/metaspace-mol-cloud/hmdb_class_of_interest.csv')
#%%

#%%
