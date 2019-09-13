#%%
#%pylab tk --no-import-all
import matplotlib
import matplotlib.pyplot as plt # For PyCharm to notice plt is defined
import json, pickle, os, gc, numpy as np, pandas as pd
from itertools import cycle, product
# Config
from analysis_pipe import filecache

PATH = '/home/lachlan/dev/notebooks/metaspace-mol-cloud'
POS, NEG = 'pos', 'neg'
# PREPROCESS:
# None - use plain cosine distances
# truncate - take N (hard-coded below) top ions based on their average cosine distance to all other mols
# normalize - stretch the 0.0001th quantile cosine distance down to 0.0001, scaling numbers below & above to match
# PREPROCESS = None
# AVG_MODE:
# 1 = naive averaging "When both molecules are present, how similar are they?"
# 2 = (incorrect) "When both molecules could potentially be present, how similar are they (filling in 0 when either is missing)?")
# 3 = "When either molecule is present, and both molecules could potentially be present, how similar are they (filling in 0 when either is missing)?"
# AVG_MODE = 3
# DATASOURCE = 'cosine' # 'deep' or 'cosine'
def get_output_filename(pol, fdr, alg, avg_mode, preprocess, func_name, ext):
    parts = [
        func_name,
        fdr and round(fdr * 100),
        alg,
        pol,
        avg_mode and {1: 'naive_avg', 2: 'broken_avg', 3: 'fixed_avg'}[avg_mode],
        preprocess and {'truncate': '_trunc', 'normalize': '_norm', None: ''}[preprocess],
    ]
    key = '_'.join(str(part) for part in parts if part)
    return f'{PATH}/charts/{key}{ext}'
FIG_WIDTH = 1920
FIG_HEIGHT = 1080

# Derived constants
matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['figure.figsize'] = (FIG_WIDTH / 100, FIG_HEIGHT / 100)

#%% Generate & save umap
if 'pivot_and_normalize_X' not in globals():
    @filecache()
    def pivot_and_normalize_X(pol, fdr, alg, avg_mode):
        raise NotImplementedError('This function is just a stub to access the cached data. Use the batch_processes2 file for putting data in the cache')

@filecache()
def run_umap(pol, fdr, alg, avg_mode, preprocess):
    import umap
    input_df = pivot_and_normalize_X(pol, fdr, alg, avg_mode)
    gc.collect()

    if preprocess == 'truncate':
        # Get top N colocalized ions
        top_ions = input_df.mean().sort_values().iloc[:5000].index
        input_df = input_df.loc[top_ions, top_ions]
    elif preprocess == 'normalize':
        base = 0.0001
        lo, hi = np.quantile(input_df, [0.0001, 1])
        assert hi - lo > base
        squeezed = np.where(input_df < lo, input_df / lo * base, (input_df - lo) / (hi - lo - base) + base)
        input_df.iloc[:, :] = np.clip(squeezed, 0, 1)
        assert isinstance(input_df, pd.DataFrame)
    elif preprocess is None:
        pass

    embedding = umap.UMAP(random_state=42, metric='precomputed').fit_transform(input_df.values)
    output_df = pd.DataFrame(embedding, columns=('x', 'y'))
    output_df['ion'] = input_df.index.values
    return output_df

# run_umap(POS, 0.05, 'cosine', 3, None)
#%% Load umap

def clip_and_join_umap(pol, fdr, alg, avg_mode, preprocess, q):
    df = run_umap(pol, fdr, alg, avg_mode, preprocess)

    xb, xlo, xhi, xt = np.quantile(df.x, [0, q, 1-q, 1])
    yb, ylo, yhi, yt = np.quantile(df.y, [0, q, 1-q, 1])
    xspan, yspan = xhi-xlo, yhi-ylo
    xmin, xmax = xlo - xspan/20, xhi + xspan/20
    ymin, ymax = ylo - yspan/20, yhi + yspan/20

    print('xmin', xmin, xmax, 'xb', xb, xt)
    print('ymin', ymin, ymax, 'yb', yb, yt)
    df = df[(df.x >= xmin) & (df.x <= xmax) & (df.y >= ymin) & (df.y <= ymax)]

    lipids_by_ion = pd.read_csv(f'{PATH}/lipids_by_ion.csv').set_index('ion')
    class_of_interest = pd.read_csv(f'{PATH}/hmdb_class_of_interest.csv').set_index('ion')
    return df.set_index('ion').join(lipids_by_ion).join(class_of_interest)

# debug = clip_and_join_umap(POS, 0.05, 'cosine', 3, None, 0.02)
#%% Plot "Is a lipid?" charts

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

def plot_lipid_chart(fdr, alg, avg_mode, preprocess, q, mpl):
    plt.figure(1).clear()
    fig, axs = plt.subplots(2,2, gridspec_kw=dict(wspace=0.025, hspace=0.25, left=0.01, top=0.95, right=0.99, bottom=0.01), num=1)
    fig.tight_layout()
    mdf_pos = clip_and_join_umap(POS, fdr, alg, avg_mode, preprocess, q)
    mdf_neg = clip_and_join_umap(NEG, fdr, alg, avg_mode, preprocess, q)

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

    for ax, df, col, title in plots:
        if mpl:
            plot_mpl(ax, df, col, title)
        else:
            plot_ds(ax, df, col, title)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])

    output = get_output_filename(None, fdr, alg, avg_mode, preprocess, 'lipids', '.png')
    print(f'Saving {output}')
    plt.savefig(f'{output}')


plot_lipid_chart(0.05, 'cosine', 3, None, 0.02, True)

#%% Plot charts highlighting individual molecule classes


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

def plot_mol_classes_chart(fdr, alg, avg_mode, preprocess, q, mpl):
    mdf_pos = clip_and_join_umap(POS, fdr, alg, avg_mode, preprocess, q)\
        .assign(class_tree=lambda df: df.super_class+'/'+df['class'],
                sub_class_tree=lambda df: df['class']+'/'+df.sub_class)
    mdf_neg = clip_and_join_umap(NEG, fdr, alg, avg_mode, preprocess, q)\
        .assign(class_tree=lambda df: df.super_class+'/'+df['class'],
                sub_class_tree=lambda df: df['class']+'/'+df.sub_class)

    plt.figure(1).clear()
    fig, axs = plt.subplots(3, 6, gridspec_kw=dict(wspace=0.025, hspace=0.25, left=0.01, top=0.95, right=0.99, bottom=0.01), num=1)
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
    for ax, is_pos, col, ind, color in zip(axs, is_poss, cols, inds, colors):
        mdf = mdf_pos if is_pos else mdf_neg
        value = mdf[col].value_counts().index[ind]
        value_str = value.replace("/", " ->\n")
        title = f'{"Positive" if is_pos else "Negative"} Mode\n{value_str}'
        if mpl:
            plot_mpl_highlight(ax, mdf, col, value, color, title)
        else:
            plot_ds_highlight(ax, mdf, col, value, color, title)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])


    output = get_output_filename(None, fdr, alg, avg_mode, preprocess, 'mol_classes', '.png')
    print(f'Saving {output}')
    plt.savefig(f'{output}')

plot_mol_classes_chart(0.05, 'cosine', 3, None, 0.02, True)

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
# ptable = ptable.apply(np.log)
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
# mdf = mdf[mdf.polarity != 'Positive mode']

def plot_mpl_rakhlin(ax, mdf, col, value, color):
    active = mdf[mdf[col] == value]
    inactive = mdf[~mdf.index.isin(active.index)]
    active_s = max(0.5, 5-len(active)**0.5/10)
    alpha = 0.75 if active_s > 3 else 0.5

    # mdf = mdf.sample(frac=1, random_state=1)
    # ax.scatter(mdf.x, mdf.y, c=[color if v == value else 'grey' for v in mdf[col]], s=[active_s if v == value else 1 for v in mdf[col]], alpha=0.5)
    ax.scatter(inactive.x, inactive.y, c=['grey'], s=1, alpha=0.1)
    ax.scatter(active.x, active.y, c=[color], s=active_s, alpha=alpha)

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
plt.savefig(f'{OUTPUT_PREFIX}_alexanders_embedding.png')

#%% Plot by hand-picked class
plt.close('all')

def plot_classes_mpl(ax, mdf, col):
    groups = mdf.groupby(col)
    categories = [('Not', 'Not', '#DDDDDD'), ('Possibly', 'Possibly', 'red'), ('Likely', 'Likely', 'green')]
    for val, label, color in categories:
        if val in groups.groups:
            grp = groups.get_group(val)
            ax.scatter(grp.x.values, grp.y.values, c=[color], s=1, label=label, alpha=0.5)
        else:
            ax.scatter([], [], c=[color], s=1, label=label, alpha=0.5)

    ax.legend(loc='upper left', markerscale=5)

for mdf, suf in [(df_pos, '_pos'), (df_neg, '_neg')]:
    plt.figure(1).clear()
    fig, axs = plt.subplots(3, 6, gridspec_kw=dict(wspace=0.025, hspace=0.25, left=0.01, top=0.95, right=0.99, bottom=0.01), num=1)
    fig.tight_layout(rect=(0,0,1,1))
    axs = axs.flatten()

    for ax, col in zip(axs, class_of_interest.columns):
        plot_classes_mpl(ax, mdf, col)
        ax.set_title(col)
        ax.set_xticks([])
        ax.set_yticks([])


    for ax in axs[len(class_of_interest.columns):]:
        ax.set_visible(False)

    print(f'Saving {OUTPUT_PREFIX}_classes{suf}.png')
    plt.savefig(f'{OUTPUT_PREFIX}_classes{suf}.png')

#%% Plot all classes
plt.close('all')
ratio_array_df = pd.read_pickle('/home/lachlan/dev/notebooks/metaspace-mol-cloud/hmdb_possibility_pivot.pickle')
mdf_pos = df_pos[['x','y']].join(ratio_array_df)
mdf_neg = df_neg[['x','y']].join(ratio_array_df)
mddf_pos = ddf_pos[['x','y']].join(ratio_array_df)
mddf_neg = ddf_neg[['x','y']].join(ratio_array_df)

categories = [('Not', 'Not', '#DDDDDD'), ('Possibly', 'Possibly', 'red'), ('Likely', 'Likely', 'green')]

for cls in ratio_array_df.columns:
    plt.figure(1).clear()
    fig, axs = plt.subplots(2, 2, gridspec_kw=dict(wspace=0.025, hspace=0.25, left=0.01, top=0.95, right=0.99, bottom=0.01), num=1)
    fig.tight_layout(rect=(0,0,1,1))
    axs = axs.flatten()

    has_data = 0
    for ax, (mdf, md) in zip(axs, [(mdf_pos, 'Cosine / Positive mode'), (mdf_neg, 'Cosine / Negative mode'),
                                   (mddf_pos, 'DL / Positive mode'), (mddf_neg, 'DL / Negative mode'),]):
        groups = mdf[['x', 'y', cls]].fillna('Not').groupby(cls)

        for val, label, color in categories:
            if val in groups.groups:
                grp = groups.get_group(val)
                ax.scatter(grp.x.values, grp.y.values, c=[color], s=1, label=label, alpha=1)
                if val != 'Not': has_data = max(has_data, len(grp))
            else:
                ax.scatter([], [], c=[color], s=1, label=label, alpha=0.5)

        ax.legend(loc='upper left', markerscale=5)
        ax.set_title(md)
        ax.set_xticks([])
        ax.set_yticks([])

    if has_data >= 20:
        print(f'Saving {OUTPUT_PREFIX}/{cls}.png')
        os.makedirs(os.path.dirname(f'{OUTPUT_PREFIX}/20_ion_threshold/{cls}.png'), exist_ok=True)
        plt.savefig(f'{OUTPUT_PREFIX}/20_ion_threshold/{cls}.png')
    # if has_data <= 5:
    #     print(f'Skipping {OUTPUT_PREFIX}/{cls}.png')



#%%