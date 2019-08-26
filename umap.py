#%%
#%pylab tk --no-import-all
import matplotlib
import matplotlib.pyplot as plt # For PyCharm to notice plt is defined
import json, pickle, os, gc, numpy as np, pandas as pd
from itertools import cycle, product
# Config
PATH = '/home/lachlan/dev/notebooks/metaspace-mol-cloud'
FDR = 0.05
FDR_PCT = round(FDR * 100)
# PREPROCESS:
# None - use plain cosine distances
# truncate - take N (hard-coded below) top ions based on their average cosine distance to all other mols
# normalize - stretch the 0.0001th quantile cosine distance down to 0.0001, scaling numbers below & above to match
PREPROCESS = None
# AVG_MODE:
# 1 = naive averaging "When both molecules are present, how similar are they?"
# 2 = (incorrect) "When both molecules could potentially be present, how similar are they (filling in 0 when either is missing)?")
# 3 = "When either molecule is present, and both molecules could potentially be present, how similar are they (filling in 0 when either is missing)?"
AVG_MODE = 3
DATASOURCE = 'cosine' # 'deep' or 'cosine'
FIG_WIDTH = 1920
FIG_HEIGHT = 1080

# Derived constants
PREPROCESS_STR = {'truncate': '_trunc', 'normalize': '_norm', None: ''}[PREPROCESS]
AVG_STR = {1: 'naive_avg', 2: 'broken_avg', 3: 'fixed_avg'}[AVG_MODE]
DATASOURCE_STR = {'cosine': ''}.get(DATASOURCE, '_' + DATASOURCE)
FN_AVG_COLOC_POS = f'{PATH}/processed_data/fdr_{FDR_PCT}/X_pos_{AVG_STR}{DATASOURCE_STR}.pickle'
FN_AVG_COLOC_NEG = f'{PATH}/processed_data/fdr_{FDR_PCT}/X_neg_{AVG_STR}{DATASOURCE_STR}.pickle'
FN_UMAP_POS = f'{PATH}/processed_data/fdr_{FDR_PCT}/umap_{DATASOURCE}{PREPROCESS_STR}_pos.pickle'
FN_UMAP_NEG = f'{PATH}/processed_data/fdr_{FDR_PCT}/umap_{DATASOURCE}{PREPROCESS_STR}_neg.pickle'
OUTPUT_PREFIX = f'{PATH}/charts/{DATASOURCE}{PREPROCESS_STR}_{FDR_PCT}'
matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['figure.figsize'] = (FIG_WIDTH / 100, FIG_HEIGHT / 100)

#%%
# Generate & save umap
def run_umap(input_fn, output_fn):
    import umap
    input_df = pd.read_pickle(input_fn)
    gc.collect()

    if PREPROCESS == 'truncate':
        # Get top N colocalized ions
        top_ions = input_df.mean().sort_values().iloc[:5000].index
        input_df = input_df.loc[top_ions, top_ions]
    elif PREPROCESS == 'normalize':
        base = 0.0001
        lo, hi = np.quantile(input_df, [0.0001, 1])
        assert hi - lo > base
        squeezed = np.where(input_df < lo, input_df / lo * base, (input_df - lo) / (hi - lo - base) + base)
        input_df.iloc[:, :] = np.clip(squeezed, 0, 1)
        assert isinstance(input_df, pd.DataFrame)

    embedding = umap.UMAP(random_state=42, metric='precomputed').fit_transform(input_df.values)
    output_df = pd.DataFrame(embedding, columns=('x', 'y'))
    output_df['ion'] = input_df.index.values
    print(f'writing {output_fn}')
    output_df.to_pickle(output_fn)


REGENERATE_UMAP = True
if REGENERATE_UMAP:
    run_umap(FN_AVG_COLOC_POS, FN_UMAP_POS)
    run_umap(FN_AVG_COLOC_NEG, FN_UMAP_NEG)
#%% Load umap

def clip_xy(df):
    ASPECT = 4/3
    xb, xlo, xhi, xt = np.quantile(df.x, [0, 0.02, 0.98, 1])
    yb, ylo, yhi, yt = np.quantile(df.y, [0, 0.02, 0.98, 1])
    xspan, yspan = xhi-xlo, yhi-ylo
    xmin, xmax = xlo - xspan/20, xhi + xspan/20
    ymin, ymax = ylo - yspan/20, yhi + yspan/20
    # overclip = xspan / yspan / ASPECT
    # print(overclip)
    # print('x', xb, xt), print(xlo, xhi), print(xmin, xmax, xspan)
    # print('y', yb, yt), print(ylo, yhi), print(ymin, ymax, yspan)
    # if overclip > 1:
    #     print('expand y')
    #     ymin -= yspan * overclip / 2 / ASPECT
    #     ymax += yspan * overclip / 2 / ASPECT
    # else:
    #     print('expand x')
    #     xmin += xspan * overclip / 2
    #     xmax -= xspan * overclip / 2
    print('xmin', xmin, xmax, 'xb', xb, xt)
    print('ymin', ymin, ymax, 'yb', yb, yt)
    # xmin, xmax = max(xmin, xb), min(xmax, xt)
    # ymin, ymax = max(ymin, yb), min(ymax, yt)
    return df[(df.x >= xmin) & (df.x <= xmax) & (df.y >= ymin) & (df.y <= ymax)]
    # return df

df_pos = clip_xy(pd.read_pickle(FN_UMAP_POS))
df_neg = clip_xy(pd.read_pickle(FN_UMAP_NEG))
lipids_by_ion = pd.read_csv(f'{PATH}/lipids_by_ion.csv').set_index('ion')
df_pos = df_pos.set_index('ion').join(lipids_by_ion)
df_neg = df_neg.set_index('ion').join(lipids_by_ion)


#%% Plot "Is a lipid?" charts
plt.close('all')
plt.figure()
# fig, axs = plt.subplots(2,2, gridspec_kw=dict(wspace=0.025, hspace=0.25, left=0.01, top=0.95, right=0.99, bottom=0.01), num=1)
fig, axs = plt.subplots(2,2, num=1)
fig.tight_layout()
mdf_pos = df_pos
mdf_neg = df_neg

plots = [
    # Plot 1 - lipid classes
    (axs[0,0], mdf_pos, 'is_lipid', 'Positive mode - Is ion a lipid?'),
    (axs[1,0], mdf_pos, 'lipid_type', 'Positive mode - Lipid type'),
    (axs[0,1], mdf_neg, 'is_lipid', 'Negative mode - Is ion a lipid?'),
    (axs[1,1], mdf_neg, 'lipid_type', 'Negative mode - Lipid type'),

    # Plot 2 -
    # (axs[0,0], mdf_pos, 'class', 'Positive mode - Molecular class'),
    # (axs[1,0], mdf_pos, 'sub_class', 'Positive mode - Molecular subclass'),
    # (axs[0,1], mdf_neg, 'class', 'Negative mode - Molecular class'),
    # (axs[1,1], mdf_neg, 'sub_class', 'Negative mode - Molecular subclass'),
]

def plot_mpl(ax, mdf, col, title):
    mdf[col] = mdf[col].astype('category')
    groups = mdf.groupby(col)
    if len(groups.groups) < 10:
        colors = plt.cm.Set1.colors
    else:
        colors = plt.cm.tab20.colors + plt.cm.tab20b.colors + plt.cm.tab20c.colors
    # categories = mdf[col].cat.categories.values
    categories = mdf[col].value_counts().index.values
    for i, (key, color) in enumerate(zip(categories, cycle(colors))):
        grp = groups.get_group(key)
        label = {'label': key} if i < 20 else {'label': ''}
        ax.scatter(grp.x.values, grp.y.values, c=[color], s=3, **label, alpha=0.75)
    # cmap = dict(zip(categories, cycle(colors)))
    print(categories)
    print(set(mdf[col].values))

    ax.legend(loc='upper left', markerscale=5)

def plot_ds(ax, mdf, col, title):
    import datashader as ds
    from datashader import transfer_functions as tf
    from datashader.colors import Sets1to3
    mdf[col] = mdf[col].astype('category')
    cvs = ds.Canvas(plot_width=960, plot_height=540)
    agg = cvs.points(mdf, 'x', 'y', ds.count_cat(col))

    color_key = dict(zip(mdf[col].cat.categories, Sets1to3))
    img = tf.shade(agg, name=title, min_alpha=64, color_key=color_key)
    img = tf.dynspread(img)
    img = img.data.view(np.uint8).reshape(*img.shape, 4)[::-1,:,:]

    ax.imshow(img)
    from matplotlib.lines import Line2D
    ax.legend([Line2D([0], [0], color=c, lw=4) for c in color_key.values()], color_key.keys(), loc='upper left', markerscale=5)

MPL = True
# MPL = False
for ax, df, col, title in plots:
    if MPL:
        plot_mpl(ax, df, col, title)
    else:
        plot_ds(ax, df, col, title)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])

print(f'Saving {OUTPUT_PREFIX}_lipids.png')
plt.savefig(f'{OUTPUT_PREFIX}_lipids.png')

#%% Plot charts highlighting individual molecule classes

mdf_pos = df_pos\
    .assign(class_tree=lambda df: df.super_class+'/'+df['class'],
            sub_class_tree=lambda df: df['class']+'/'+df.sub_class)
mdf_neg = df_neg\
    .assign(class_tree=lambda df: df.super_class+'/'+df['class'],
            sub_class_tree=lambda df: df['class']+'/'+df.sub_class)

def plot_mpl_highlight(ax, mdf, col, value, color, title):
    # categories = mdf[col].cat.categories.values
    active = mdf[mdf[col] == value]
    inactive = mdf[~mdf.index.isin(active.index)]
    active_s = max(0.5, 5-len(active)**0.5/10)
    alpha = 0.75 if active_s > 3 else 0.5

    ax.scatter(inactive.x, inactive.y, c=['grey'], s=1, alpha=0.1, label='Other')
    ax.scatter(active.x, active.y, c=[color], s=active_s, alpha=alpha, label=value)

def plot_ds_highlight(ax, mdf, col, value, color, title):
    global img
    import datashader as ds
    from datashader import transfer_functions as tf
    mdf = mdf[['x','y']].assign(highlight=(mdf[col]==value).astype('category'))
    cvs = ds.Canvas(plot_width=410, plot_height=260)
    agg = cvs.points(mdf, 'x', 'y', ds.count_cat('highlight'))

    img = tf.shade(agg, name=title, min_alpha=64, color_key={False: 'grey', True: tuple(np.array(color)*255)})
    img = tf.dynspread(img)
    img = img.data.view(np.uint8).reshape(*img.shape, 4)[::-1, :, :]

    ax.imshow(img)

plt.figure(1).clear()
fig, axs = plt.subplots(4, 6, gridspec_kw=dict(wspace=0.025, hspace=0.25, left=0.01, top=0.95, right=0.99, bottom=0.01), num=1)
fig.tight_layout(rect=(0,0,1,1))
axs = [*axs[:, :3].flatten(), *axs[:, 3:].flatten()]
is_poss = [True] * 12 + [False] * 12
# cols = (['class_tree'] * 9 + ['sub_class_tree'] * 3) * 2
# cols = (['super_class'] * 6 + ['class_tree'] * 6) * 2
# inds = ([*range(9), *range(3)]) * 2
cols = (['class_tree'] * 12) * 2
inds = ([*range(12)]) * 2
# colors = ([*(plt.cm.Set1.colors + plt.cm.Set2.colors + plt.cm.Set3.colors)][:12]) * 2
colors = cycle(plt.cm.tab10.colors[:4] + plt.cm.tab10.colors[5:7]) # Only non-greyish colors
MPL = True
# MPL = False
for ax, is_pos, col, ind, color in zip(axs, is_poss, cols, inds, colors):
    mdf = mdf_pos if is_pos else mdf_neg
    value = mdf[col].value_counts().index[ind]
    value_str = value.replace("/", " ->\n")
    title = f'{"Positive" if is_pos else "Negative"} Mode\n{value_str}'
    if MPL:
        plot_mpl_highlight(ax, mdf, col, value, color, title)
    else:
        plot_ds_highlight(ax, mdf, col, value, color, title)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])


print(f'Saving {OUTPUT_PREFIX}_mol_classes.png')
plt.savefig(f'{OUTPUT_PREFIX}_mol_classes.png')

#%% Bespoke plot

plt.close('all')
# fig, axs = plt.subplots(1, 2, gridspec_kw=dict(wspace=0.025, hspace=0.025, left=0.01, top=0.95, right=0.99, bottom=0.01), num=1)
fig, ax = plt.figure(), plt.gca()
fig.tight_layout(rect=(0,0,1,1))

def plot_mpl_bespoke(ax, mdf, cats, explicit_colors):
    from datashader.colors import Sets1to3
    color_key = dict(zip(cats, Sets1to3))
    color_key.update(explicit_colors)
    s = max(0.5, 20 - len(mdf) ** 0.5 / 10)
    ax.scatter(mdf.x, mdf.y, c=[color_key[cat] for cat in mdf.cat], s=s, alpha=0.75, marker='.')

    from matplotlib.lines import Line2D
    ax.legend([Line2D([0], [0], color=color_key[c], lw=4) for c in cats], cats, loc='upper left', markerscale=5)


def plot_ds_bespoke(ax, mdf, cats, explicit_colors):
    import datashader as ds
    from datashader import transfer_functions as tf
    from datashader.colors import Sets1to3
    cvs = ds.Canvas(plot_width=960, plot_height=640)
    agg = cvs.points(mdf, 'x', 'y', ds.count_cat('cat'))

    color_key = dict(zip(cats, Sets1to3))
    color_key.update(explicit_colors)
    img = tf.shade(agg, min_alpha=128, color_key=color_key)
    img = tf.dynspread(img, threshold=0.5, max_px=1000)
    img = img.data.view(np.uint8).reshape(*img.shape, 4)[::-1,:,:]

    ax.imshow(img)
    from matplotlib.lines import Line2D
    ax.legend([Line2D([0], [0], color=color_key[c], lw=4) for c in cats], cats, loc='upper left', markerscale=5)

cats = ['Glycerophospholipids', 'Glycerolipids', 'Fatty Acyls', 'Flavonoids']

def get_cat_idx(row): return row['class'] if row['class'] in cats else 'Other metabolites'
mdf = df_pos[['x', 'y']].assign(cat=df_pos.apply(get_cat_idx, axis=1).astype('category'))
MPL = True
# MPL = False
if MPL:
    plot_mpl_bespoke(ax, mdf, [*cats, 'Other metabolites'], {'Other metabolites': 'lightgray'})
else:
    plot_ds_bespoke(ax, mdf, [*cats, 'Other metabolites'], {'Other metabolites': 'lightgray'})
ax.set_xticks([])
ax.set_yticks([])

print(f'Saving {OUTPUT_PREFIX}_bespoke.png')
plt.savefig(f'{OUTPUT_PREFIX}_bespoke.png')

#%% Alexander Rakhlin's data

df = pd.read_csv(f'{PATH}/raw_data/preds_pi_for_lachlan_fdr5.csv')
df['source'] = df['sourceSf'] + df['sourceAdduct']
df['target'] = df['targetSf'] + df['targetAdduct']
df = df.drop(labels=['sourceSf', 'sourceAdduct', 'targetSf', 'targetAdduct'], axis=1)
df[['source', 'target']] = np.sort(df[['source', 'target']].values, axis=1)
groups = df.groupby(['source', 'target'])
df2 = groups.aggregate(np.median).reset_index()
counts = pd.concat([df2['source'], df2['target']]).value_counts()
freq = 1700
print(f'# of frequent ions. those that appear >{freq} times in this case: {sum(counts > freq)}')
basis = list(counts[counts > freq].index)
idx = df2['source'].isin(basis) | df2['target'].isin(basis)
data = df2[idx]

# make sure that all target ions belong to basis - rearrange when needed
ixf = ~data['target'].isin(basis)
data.loc[ixf, ['source', 'target']] = data[['source', 'target']].values[ixf, ::-1]
assert sum(~data['target'].isin(basis)) == 0
# create distance features
ptable = pd.pivot_table(data, index='source', columns='target', values='score')
print(ptable.shape)

assert sum(ptable.columns.isin(basis)) == len(basis) == ptable.shape[1]
avnc = int(ptable.notna().sum(1).mean()) # average non-sparse coordinates per row
print(f'distance table in reduced basis: shape {ptable.shape}, average non-sparse coordinates per row {avnc}, non-sparse ratio {int(avnc / ptable.shape[1] * 100)}%')

btbl = pd.concat([df2.groupby('source').apply(len), df2.groupby('target').apply(len)]).groupby(level=0).apply(sum)
shape = (btbl.shape[0], btbl.shape[0])
print(f'full distance table: shape {shape}, average non-sparse coordinates per row {int(btbl.mean())}, non-sparse ratio {int(btbl.mean() / shape[0] * 100)}%')
## Embedding
ptable = ptable.apply(np.log)
coo = ptable.to_sparse().to_coo()
# U-map embedding
MIN_DIST = 1.0
N_NEIGHBORS = 600
N_COMPONENTS = 2
METRIC = 'manhattan'
RANDOM_STATE = 0

import umap
embd = umap.UMAP(metric=METRIC, n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST,
                 n_components=N_COMPONENTS, random_state=RANDOM_STATE).fit_transform(coo)
#%%

import re
mdf = (pd.DataFrame(embd, index=ptable.index, columns=['x','y'], copy=True)
       .join(lipids_by_ion)
       .assign(polarity=lambda df: ['Positive mode' if re.search('\+(H|Na|K)$', s) else 'Negative mode' for s in ptable.index],
               class_tree=lambda df: df.super_class+'/'+df['class'],
               sub_class_tree=lambda df: df['class']+'/'+df.sub_class))
# mdf = mdf[mdf.y > -20]
# mdf = mdf[mdf.polarity == 'Positive mode']

def plot_mpl_rakhlin(ax, mdf, col, value, color):
    active = mdf[mdf[col] == value]
    inactive = mdf[~mdf.index.isin(active.index)]
    active_s = max(0.5, 5-len(active)**0.5/10)
    alpha = 0.75 if active_s > 3 else 0.5

    # mdf = mdf.sample(frac=1, random_state=1)
    # ax.scatter(mdf.x, mdf.y, c=[color if v == value else 'grey' for v in mdf[col]], s=[active_s if v == value else 1 for v in mdf[col]], alpha=0.5)
    ax.scatter(-inactive.x, -inactive.y, c=['grey'], s=1, alpha=0.1)
    ax.scatter(-active.x, -active.y, c=[color], s=active_s, alpha=alpha)

plt.figure(1).clear()
fig, axs = plt.subplots(3, 4, gridspec_kw=dict(wspace=0.025, hspace=0.25, left=0.01, top=0.95, right=0.99, bottom=0.01), num=1)
fig.tight_layout(rect=(0,0,1,1))
axs = axs.flatten()
is_poss = [True] * 12
cols = ['polarity'] + ['class_tree'] * 12
inds = range(12)
colors = cycle(plt.cm.tab10.colors[:4] + plt.cm.tab10.colors[5:7]) # Only non-greyish colors
for ax, is_pos, col, ind, color in zip(axs, is_poss, cols, inds, colors):
    value = mdf[col].value_counts().index[ind]
    value_str = value.replace("/", " ->\n")
    plot_mpl_rakhlin(ax, mdf, col, value, color)
    ax.set_title(value_str)
    ax.set_xticks([])
    ax.set_yticks([])


print(f'Saving {OUTPUT_PREFIX}_alexanders_embedding.png')
plt.savefig(f'{OUTPUT_PREFIX}_alexanders_embedding_outlier.png')

#%%