#%%
import os; os.chdir('/home/lachlan/dev/notebooks/metaspace-mol-cloud/specific_classes')
import matplotlib, matplotlib.colors
# %matplotlib tk
import matplotlib.pyplot as plt
import os, numpy as np, pandas as pd
# Config
PATH = '.'
OUTPUT_PATH = f'{PATH}/output'

FIG_WIDTH = 1600
FIG_HEIGHT = 850

matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['figure.figsize'] = (FIG_WIDTH / 100, FIG_HEIGHT / 100)

#%% Load umap
def clip_xy(df, bottom=0, top=0, left=0, right=0):
    """Remove outlying N% of points from top, left, bottom and right of graph"""
    xb, xlo, xhi, xt = np.percentile(df.x, [0, left, 100-right, 100])
    yb, ylo, yhi, yt = np.percentile(df.y, [0, bottom, 100-top, 100])
    xspan, yspan = xhi-xlo, yhi-ylo
    xmin, xmax = xlo - xspan/20, xhi + xspan/20
    ymin, ymax = ylo - yspan/20, yhi + yspan/20
    # print('xmin', xmin, xmax, 'xb', xb, xt)
    # print('ymin', ymin, ymax, 'yb', yb, yt)
    return df[(df.x >= xmin) & (df.x <= xmax) & (df.y >= ymin) & (df.y <= ymax)]


ratio_array_df = pd.read_pickle(f'{PATH}/hmdb_possibility_pivot.pickle')
cosine_pos = clip_xy(pd.read_pickle(f'{PATH}/umap_cosine_pos.pickle'), top=1, bottom=1).set_index('ion').join(ratio_array_df)
cosine_neg = clip_xy(pd.read_pickle(f'{PATH}/umap_cosine_neg.pickle'), right=1).set_index('ion').join(ratio_array_df)
deep_pos = clip_xy(pd.read_pickle(f'{PATH}/umap_deep_pos.pickle'), top=2, right=2).set_index('ion').join(ratio_array_df)
deep_neg = clip_xy(pd.read_pickle(f'{PATH}/umap_deep_neg.pickle'), left=2).set_index('ion').join(ratio_array_df)

#%% Load umap

# @filecache(lambda df, cls, ds, r: f'{cls[cls.rindex("/")+1:]}_{ds}_{r}')
def get_n_neighbours(df, cls, ds, r):
    from scipy.spatial.ckdtree import cKDTree
    global points, tree
    points = np.array([df.x, df.y]).T
    points -= np.min(points, axis=0, keepdims=True)
    points /= np.max(points, axis=0, keepdims=True)
    tree = cKDTree(points)
    return pd.Series([*map(len, tree.query_ball_point(points, r))], index=df.index)



#%% Plot specific classes
# plt.close('all')

# classes_to_plot = ratio_array_df.columns
classes_to_plot = [
    'Lipids and lipid-like molecules/Glycerolipids',
    'Lipids and lipid-like molecules/Glycerolipids/Triradylcglycerols',
    'Lipids and lipid-like molecules/Glycerophospholipids',
    'Lipids and lipid-like molecules/Glycerophospholipids/Glycerophosphoethanolamines',
]

categories = ['Not', 'Possibly', 'Likely']

def get_style(mol_class, datasource, category, df):
    category_marker_styles = dict([
        ('Not', dict(s=32, alpha=0.7, c=['#888888'], marker='o', linewidths=0)),
        ('Possibly', dict(s=128, alpha=0.7, c=['tab:red'], marker='o', linewidths=0)),
        ('Likely', dict(s=128, alpha=0.7, c=['tab:green'], marker='o', linewidths=0)),
    ])
    style = category_marker_styles[category]
    print(len(df))
    if datasource == 'Cosine / Positive mode' and category == 'Not':
        style['alpha'] = 0.2
    elif datasource == 'Cosine / Negative mode' and category == 'Not':
        style['alpha'] = 0.4
    elif len(df[df[mol_class] != 'Not']) > 500:
        style['alpha'] = 0.5

    return style

def plot(cls, mdf, md):
    plt.figure(1).clear()
    fig, ax = plt.subplots(1, 1, gridspec_kw=dict(wspace=0.025, hspace=0.06, left=0.001, top=0.999, right=0.999, bottom=0.001), num=1)
    fig.tight_layout(rect=(0,0,1,1))

    mdf = mdf[['x', 'y', cls]].fillna('Not')
    groups = mdf.groupby(cls)

    for val in categories:
        if val in groups.groups:
            grp = groups.get_group(val)
            style = get_style(cls, md, val, mdf)
            ax.scatter(grp.x.values, grp.y.values, **style)

    # ax.legend(loc='upper left', markerscale=5)
    # ax.set_title(md)
    ax.set_xticks([])
    ax.set_yticks([])
    name = md.replace(' / ', ' ') + ' ' + (cls[cls.rindex('/') + 1:] if '/' in cls else cls)
    print(f'Saving {OUTPUT_PATH}/{name}.png')
    os.makedirs(os.path.dirname(f'{OUTPUT_PATH}/{name}.png'), exist_ok=True)
    plt.savefig(f'{OUTPUT_PATH}/{name}.png')

plots = [
    ('Lipids and lipid-like molecules/Glycerolipids', cosine_pos, 'Cosine / Positive mode'),
    ('Lipids and lipid-like molecules/Glycerolipids/Triradylcglycerols', cosine_pos, 'Cosine / Positive mode'),
    ('Lipids and lipid-like molecules/Glycerolipids', cosine_neg, 'Cosine / Negative mode'),
    ('Lipids and lipid-like molecules/Glycerolipids', deep_pos, 'DL / Positive mode'),
    ('Lipids and lipid-like molecules/Glycerophospholipids', cosine_pos, 'Cosine / Positive mode'),
    ('Lipids and lipid-like molecules/Glycerophospholipids/Glycerophosphoethanolamines', cosine_pos, 'Cosine / Positive mode'),
    ('Lipids and lipid-like molecules/Glycerophospholipids', cosine_neg, 'Cosine / Negative mode'),
    ('Lipids and lipid-like molecules/Glycerophospholipids', deep_pos, 'DL / Positive mode'),
]

for cls, mdf, md in plots:
    # for (mdf, md) in [(cosine_pos, 'Cosine / Positive mode'), (cosine_neg, 'Cosine / Negative mode'),
    #                                (deep_pos, 'DL / Positive mode'), (deep_neg, 'DL / Negative mode'),]:
    plot(cls, mdf, md)

# plot(classes_to_plot[0], cosine_pos, 'Cosine / Positive mode')

# %%