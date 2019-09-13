
#%% Output for cytoscape
sources = [
    # ('all', merged_coloc.groupby(['source','target'])),
    # ('above_0', merged_coloc[lambda df: df.cosine > 0].groupby(['source','target'])),
    # ('above_50', merged_coloc[lambda df: df.cosine > 0.5].groupby(['source','target'])),
    # ('above_70', merged_coloc[lambda df: df.cosine > 0.7].groupby(['source','target']))
    ('above_90', merged_coloc[lambda df: df.cosine > 0.9].groupby(['source', 'target']))
]
dfs = []
for prefix, groupby in sources:
    dfs.append(groupby.cosine.mean().rename(f'{prefix}_cosine_mean'))
    dfs.append(groupby.cosine.sum().rename(f'{prefix}_cosine_sum'))
    dfs.append(groupby.cosine.count().rename(f'{prefix}_cosine_count'))

mega_table = pd.concat(dfs, axis=1).reset_index()
mega_table.to_csv('/home/lachlan/dev/notebooks/metaspace-mol-cloud/flat_coloc_fdr_10_above_90pct.csv', index=False)
