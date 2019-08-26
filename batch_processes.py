#%%
import json, pickle, os, gc
from datetime import datetime

import numpy as np
import pandas as pd
from metaspace.sm_annotation_utils import SMInstance, SMDataset
from concurrent.futures import ProcessPoolExecutor

PATH = '/home/lachlan/dev/notebooks/metaspace-mol-cloud'
FDR = 0.05
FDR_PCT = round(FDR * 100)
# AVG_MODE:
# 1 = naive averaging "When both molecules are present, how similar are they?"
# 2 = (incorrect) "When both molecules could potentially be present, how similar are they (filling in 0 when either is missing)?")
# 3 = "When either molecule is present, and both molecules could potentially be present, how similar are they (filling in 0 when either is missing)?"
AVG_MODE = 3
AVG_STR = {1: 'naive_avg', 2: 'broken_avg', 3: 'fixed_avg'}[AVG_MODE]
DATASOURCE = 'deep' # 'deep' or 'cosine'
DATASOURCE_STR = {'deep':'_deep', 'cosine': ''}[DATASOURCE]
FN_MS_DATA_POS = f'{PATH}/raw_data/metaspace_data_pos_fdr_{FDR_PCT}.pickle'
FN_MS_DATA_NEG = f'{PATH}/raw_data/metaspace_data_neg_fdr_{FDR_PCT}.pickle'
FN_RAW_COLOC_POS = f'{PATH}/raw_data/coloc_pos_hmdb{DATASOURCE_STR}_fdr_{FDR_PCT}.pickle'
FN_RAW_COLOC_NEG = f'{PATH}/raw_data/coloc_neg_hmdb{DATASOURCE_STR}_fdr_{FDR_PCT}.pickle'
FN_DENOM_POS = f'{PATH}/raw_data/denom_{AVG_STR}{DATASOURCE_STR}_pos_fdr_{FDR_PCT}.pickle'
FN_DENOM_NEG = f'{PATH}/raw_data/denom_{AVG_STR}{DATASOURCE_STR}_neg_fdr_{FDR_PCT}.pickle'
FN_AVG_COLOC_POS = f'{PATH}/processed_data/fdr_{FDR_PCT}/X_pos_{AVG_STR}{DATASOURCE_STR}.pickle'
FN_AVG_COLOC_NEG = f'{PATH}/processed_data/fdr_{FDR_PCT}/X_neg_{AVG_STR}{DATASOURCE_STR}.pickle'

#%%

#%%
def get_coloc_matrix(anns):
    for retry in range(3):
        try:
            from sm.engine.png_generator import ImageStoreServiceWrapper
            from sklearn.metrics.pairwise import pairwise_kernels
            from scipy.ndimage import median_filter

            anns = sorted(anns, key=lambda ann: -ann['msmScore'])[:2000]
            cnt = len(anns)
            img_ids = [ann['isotopeImages'][0]['url'][-32:] for ann in anns]
            img_names = [(ann['sumFormula'] + ann['adduct']) for ann in anns]
            img_svc = ImageStoreServiceWrapper('https://metaspace2020.eu/')
            images, mask, (h, w) = img_svc.get_ion_images_for_analysis('fs', img_ids, max_mem_mb=2048, hotspot_percentile=100)
            # Discard everything below 50% intensity
            images[images < np.quantile(images, 0.5, axis=1, keepdims=True)] = 0
            filtered_images = median_filter(images.reshape((cnt, h, w)), (1, 3, 3)).reshape((cnt, h * w))

            distance_matrix = pairwise_kernels(filtered_images, metric='cosine')
            df = pd.DataFrame(distance_matrix, index=img_names, columns=img_names, dtype=np.float32)
            df.rename_axis(index='source', columns='target', inplace=True)
            df.reset_index(inplace=True)
            unpivoted = df.melt(id_vars='source', value_name='cosine')
            return unpivoted[lambda df: (df.source < df.target) & (df.cosine > 0.01)]  # Reduce to uni-directional links & do initial thresholding
        except Exception as ex:
            print(ex)
    return None

def get_ds_data(ds_id):
    global FDR

    sm = SMInstance()
    anns = sm._gqclient.getAnnotations({'database': 'HMDB-v4', 'fdrLevel': FDR,
                                        'hasNeutralLoss': False, 'hasChemMod': False, 'hasHiddenAdduct': False},
                                       {'ids': ds_id})
    if len(anns) > 2:
        coloc = get_coloc_matrix(anns)
        mzs = np.array([ann['mz'] for ann in anns])
        mz_range = np.min(mzs), np.max(mzs)
        mz_dict = dict([(ann['sumFormula'] + ann['adduct'], ann['mz']) for ann in anns])
        return coloc, mz_range, mz_dict
    return None

#%% Process, Concat, Save colocalizations
sm = SMInstance()
# sm.login(**json.loads(open('/home/lachlan/.metaspace.json')))

datasets = sm.datasets()
def fetch_data_from_metaspace(is_pos, coloc_filename, data_filename):
    all_coloc = []
    all_ranges = []
    all_mz_dict = {}
    ion_present_in_ds = set()
    datasets = [SMDataset(info, sm._gqclient) for info in sm._gqclient.getDatasets({'polarity': 'POSITIVE' if is_pos else 'NEGATIVE'})]
    ds_ids = []
    with ProcessPoolExecutor(8) as ex:
        for i, result in enumerate(ex.map(get_ds_data, [ds.id for ds in datasets])):
            if i % 100 == 0: print("pos" if is_pos else "neg", i, 'of', len(datasets))
            if result is not None:
                ds_ids.append(datasets[i].id)
                coloc, mz_range, mz_dict = result
                n = len(all_coloc) # Keep a consistent index of which non-empty DS we are at
                all_coloc.append(coloc)
                all_ranges.append(mz_range)
                all_mz_dict.update(mz_dict)
                ion_present_in_ds.update((n, mol) for mol in coloc.source)
                ion_present_in_ds.update((n, mol) for mol in coloc.target)

    all_coloc = pd.concat(all_coloc, ignore_index=True)
    print('saving', coloc_filename)
    all_coloc.to_pickle(coloc_filename)
    del all_coloc; gc.collect()
    pickle.dump((ds_ids, all_ranges, all_mz_dict, ion_present_in_ds), open(data_filename, 'wb'))

if DATASOURCE == 'cosine':
    fetch_data_from_metaspace(True, FN_RAW_COLOC_POS, FN_MS_DATA_POS)
    fetch_data_from_metaspace(False, FN_RAW_COLOC_NEG, FN_MS_DATA_NEG)
#%% Load non-colocalization data
if DATASOURCE == 'cosine':
    pos_ds_ids, pos_ranges, pos_mz_dict, pos_ion_present_in_ds = pickle.load(open(FN_MS_DATA_POS, 'rb'))
    neg_ds_ids, neg_ranges, neg_mz_dict, neg_ion_present_in_ds = pickle.load(open(FN_MS_DATA_NEG, 'rb'))
#%% Filter non-colocalization data to match external data source
def filter_ext_data(ext_coloc, ms_data, output_fn):
    global i, ext_ds_ids, common_ds_idxs, ds_ids, mz_ranges
    ds_ids, mz_ranges, mz_dict, ion_present_in_ds = pickle.load(open(ms_data, 'rb'))

    coloc = ext_coloc[ext_coloc.datasetId.isin(ds_ids)][['source','target','score']]
    ext_ds_ids_set = set(ext_coloc.datasetId.unique())
    common_ds_idxs = [i for i, id in enumerate(ds_ids) if id in ext_ds_ids_set]
    common_ds_idxs_set = set(common_ds_idxs)
    ext_ds_ids = [ds_ids[i] for i in common_ds_idxs]
    ext_mz_ranges = [mz_ranges[i] for i in common_ds_idxs]
    ext_ion_present_in_ds = [(i, ion) for i, ion in ion_present_in_ds if i in common_ds_idxs_set]
    fixed_coloc = coloc.assign(score=1 - coloc.score / np.nanmax(coloc.score))

    fixed_coloc.to_pickle(output_fn)
    return ext_ds_ids, ext_mz_ranges, mz_dict, ext_ion_present_in_ds

if DATASOURCE == 'deep':
    ext_coloc = (pd.read_csv(f'{PATH}/raw_data/preds_pi_for_lachlan_fdr5.csv')
                 .assign(source=lambda df: df.sourceSf + df.sourceAdduct,
                         target=lambda df: df.targetSf + df.targetAdduct))
    pos_ds_ids, pos_ranges, pos_mz_dict, pos_ion_present_in_ds = filter_ext_data(ext_coloc, FN_MS_DATA_POS, FN_RAW_COLOC_POS)
    neg_ds_ids, neg_ranges, neg_mz_dict, neg_ion_present_in_ds = filter_ext_data(ext_coloc, FN_MS_DATA_NEG, FN_RAW_COLOC_NEG)
#%% Make table of how many potential colocalizations actually can exist (This is the denominator of the averaging step)
def make_mz_df_index(mz_ranges, mz_dict):
    mz_df = (pd.DataFrame(mz_dict.items(), columns=['ion','mz']))
    mz_df['ds_idxs'] = mz_df.mz.apply(lambda mz: np.packbits((mz_ranges[:, 0] <= mz) & (mz_ranges[:, 1] >= mz)))
    mz_df.set_index('ion', inplace=True)
    return mz_df

def make_present_index(ion_present_in_ds, ds_cnt, ion_idx):
    global present_idx
    present_idx = (pd.DataFrame(ion_present_in_ds, columns=['ds','ion'])
                   .assign(cnt=True)
                   .pivot_table(values='cnt', index='ds', columns='ion', fill_value=False)
                   .reindex(index=list(range(ds_cnt)), fill_value=False))
    table = np.zeros((len(ion_idx), int(np.ceil(ds_cnt/8))), dtype=np.uint8)
    for i, ion in enumerate(ion_idx):
        table[i, :] = np.packbits(list(present_idx[ion].values)) if ion in present_idx.columns else 0
    return table


def make_mz_coloc_df(mz_df, ion_present_table):
    import numba
    POPCOUNT_TABLE8 = np.array([bin(x).count("1") for x in range(256)], dtype=int)
    mz_coloc_df = pd.DataFrame(0, index=mz_df.index, columns=mz_df.index, dtype=np.uint16)

    @numba.jit(nopython=True, parallel=True, nogil=True)
    def do_rows(_mz_coloc_df, _ds_idxs, _ion_present_table, _avg_mode, start, end):
        for i in numba.prange(start, end):
            _mz_coloc_df[i, i] = np.sum(POPCOUNT_TABLE8[_ds_idxs[i]])
            for j in range(min(i, len(_ds_idxs))):
                mask = np.bitwise_and(_ds_idxs[i], _ds_idxs[j])
                if _avg_mode == 3:
                    either_ion_is_present = np.bitwise_or(_ion_present_table[i], _ion_present_table[j])
                    mask = np.bitwise_and(mask, either_ion_is_present)

                val = np.sum(POPCOUNT_TABLE8[mask])
                _mz_coloc_df[i, j] = _mz_coloc_df[j, i] = val

    for i in range(0, len(mz_df), 1000):
        i_end = min(i + 1000, len(mz_df))
        print(datetime.now().ctime(), i, 'to', i_end, 'out of', len(mz_df))
        do_rows(mz_coloc_df.values, list(mz_df.ds_idxs.values), ion_present_table, AVG_MODE, i, i_end)
    return mz_coloc_df

def generate_denom_matrix(mz_ranges, mz_dict, ion_present_in_ds, filename):
    global mz_df, ion_present_table, mz_coloc_df
    print(f'generate_denom_matrix{repr((len(mz_ranges), len(mz_dict), len(ion_present_in_ds), filename))}')
    mz_df = make_mz_df_index(mz_ranges, mz_dict)
    ion_present_table = make_present_index(ion_present_in_ds, len(mz_ranges), mz_df.index)
    mz_coloc_df = make_mz_coloc_df(mz_df, ion_present_table)
    mz_coloc_df.to_pickle(filename)


if AVG_MODE in (2, 3):
    # Molecules' m/z values can vary by a small amount based on resolving power :( Widen the m/z window slightly to catch molecules that lie near the border
    fixed_pos_ranges = np.array(pos_ranges) + [[-0.0025, 0.0025]]
    fixed_neg_ranges = np.array(neg_ranges) + [[-0.0025, 0.0025]]
    generate_denom_matrix(fixed_pos_ranges, pos_mz_dict, pos_ion_present_in_ds, FN_DENOM_POS)
    generate_denom_matrix(fixed_neg_ranges, neg_mz_dict, neg_ion_present_in_ds, FN_DENOM_NEG)

#%% Finally make the comparison matrix
def pivot_and_normalize_X(raw_coloc_filename, denom_filename, output_filename):
    global pre_avg, fit_denom
    print(f'pivot_and_normalize_X{repr((raw_coloc_filename, denom_filename, output_filename))}')
    gc.collect()
    aggfunc = 'sum' if AVG_MODE in (2,3) else 'mean'
    X = pd.read_pickle(raw_coloc_filename)
    value_col = 'score' if 'score' in X.columns else 'cosine'
    X = pd.concat([X, X.rename(columns={'source': 'target', 'target': 'source'})], ignore_index=True, sort=True) \
        .pivot_table(index='source', columns='target', values=value_col, aggfunc=aggfunc)
    gc.collect()
    if AVG_MODE in (2,3):
        denom = pd.read_pickle(denom_filename)
        pre_avg = X
        fit_denom = denom.reindex_like(X)
        fit_denom[fit_denom == 0] = 1
        X = X / denom.reindex_like(X)
        X.values[X.values > 1] = 1

    X.fillna(0, inplace=True)
    # Convert sum-similarity to something resembling a distance
    X = 1 - X
    for i in range(len(X)): X.iloc[i, i] = 0
    print('writing', output_filename)
    X.to_pickle(output_filename)
    del X; gc.collect()

pivot_and_normalize_X(FN_RAW_COLOC_POS, FN_DENOM_POS, FN_AVG_COLOC_POS)
pivot_and_normalize_X(FN_RAW_COLOC_NEG, FN_DENOM_NEG, FN_AVG_COLOC_NEG)

#%% Make apples-to-apples comparison by only including ions present in both
def make_intersection_df(avg_coloc_cosine_fn, avg_coloc_deep_fn, cosine_output_fn, deep_output_fn):
    cosine_coloc = pd.read_pickle(avg_coloc_cosine_fn)
    deep_coloc = pd.read_pickle(avg_coloc_deep_fn)
    common_ions = set(cosine_coloc.index).intersection(deep_coloc.index)
    cos_count = np.count_nonzero(cosine_coloc != 1) - len(cosine_coloc)
    dee_count = np.count_nonzero(deep_coloc != 1) - len(deep_coloc)

    cosine_coloc = cosine_coloc.loc[common_ions, common_ions]
    deep_coloc = deep_coloc.loc[common_ions, common_ions]

    clipped_cos_count = np.count_nonzero(cosine_coloc != 1) - len(cosine_coloc)
    clipped_dee_count = np.count_nonzero(deep_coloc != 1) - len(deep_coloc)

    uncommon_mask = ~((cosine_coloc != 1) & (deep_coloc != 1))
    cosine_coloc.values[uncommon_mask] = 1
    deep_coloc.values[uncommon_mask] = 1

    filtered_cos_count = np.count_nonzero(cosine_coloc != 1) - len(cosine_coloc)
    filtered_dee_count = np.count_nonzero(deep_coloc != 1) - len(deep_coloc)

    print('\tInitial count -> common ion count -> common pair count')
    print(f'cosine\t{cos_count:8}\t{clipped_cos_count:8}\t{filtered_cos_count:8}')
    print(f'deep  \t{dee_count:8}\t{clipped_dee_count:8}\t{filtered_dee_count:8}')

    cosine_coloc.to_pickle(cosine_output_fn)
    deep_coloc.to_pickle(deep_output_fn)

make_intersection_df(f'{PATH}/processed_data/fdr_20/X_pos_{AVG_STR}.pickle',
               f'{PATH}/processed_data/fdr_5/X_pos_{AVG_STR}_deep.pickle',
               f'{PATH}/processed_data/fdr_20/X_pos_{AVG_STR}_cosine_mutual.pickle',
               f'{PATH}/processed_data/fdr_20/X_pos_{AVG_STR}_deep_mutual.pickle')

make_intersection_df(f'{PATH}/processed_data/fdr_20/X_neg_{AVG_STR}.pickle',
               f'{PATH}/processed_data/fdr_5/X_neg_{AVG_STR}_deep.pickle',
               f'{PATH}/processed_data/fdr_20/X_neg_{AVG_STR}_cosine_mutual.pickle',
               f'{PATH}/processed_data/fdr_20/X_neg_{AVG_STR}_deep_mutual.pickle')


#%% Expand DFs to the same shape
def make_expanded_df(avg_coloc_cosine_fn, avg_coloc_deep_fn, cosine_output_fn, deep_output_fn):
    cosine_coloc = pd.read_pickle(avg_coloc_cosine_fn)
    deep_coloc = pd.read_pickle(avg_coloc_deep_fn)
    from collections import OrderedDict
    from itertools import chain
    all_ions = pd.Index(OrderedDict.fromkeys(chain(cosine_coloc.index, deep_coloc.index)))

    cosine_coloc = cosine_coloc.reindex(index=all_ions, columns=all_ions, fill_value=1)
    deep_coloc = deep_coloc.reindex(index=all_ions, columns=all_ions, fill_value=1)

    cosine_coloc.values[np.diag_indices(len(all_ions))] = 0
    deep_coloc.values[np.diag_indices(len(all_ions))] = 0

    cosine_coloc.to_pickle(cosine_output_fn)
    deep_coloc.to_pickle(deep_output_fn)

make_expanded_df(f'{PATH}/processed_data/fdr_5/X_pos_{AVG_STR}.pickle',
               f'{PATH}/processed_data/fdr_5/X_pos_{AVG_STR}_deep.pickle',
               f'{PATH}/processed_data/fdr_5/X_pos_{AVG_STR}_cosine_expand.pickle',
               f'{PATH}/processed_data/fdr_5/X_pos_{AVG_STR}_deep_expand.pickle')

make_expanded_df(f'{PATH}/processed_data/fdr_5/X_neg_{AVG_STR}.pickle',
               f'{PATH}/processed_data/fdr_5/X_neg_{AVG_STR}_deep.pickle',
               f'{PATH}/processed_data/fdr_5/X_neg_{AVG_STR}_cosine_expand.pickle',
               f'{PATH}/processed_data/fdr_5/X_neg_{AVG_STR}_deep_expand.pickle')
# %%

