#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as ss
import seaborn as sns
import os
import math
from data.sequencing import data_import
from data.sequencing import data_analysis

# Plot settings
plt.rcParams.update({'font.size': 15, 'svg.fonttype' : 'none', 
                     'figure.figsize' : [15,12], 'xtick.major.size' : 10, 
                     'ytick.major.size' : 10, 'lines.markersize' : 12})


# In[3]:


queries = {'cat' : 'ACCGGTAAGGCACAGAAACG',
           'dna' : 'ATGTGTAAGGCACAAAAACG',
           'frame' : 'ATACGCAAGGAACAAAAACG',
           'building' : 'ATTTGCAAGGAACAAAAACG',
          'lego' : 'ACCTGTAAGGCACAGAAACG'}

bls = ['bl1', 'bl2', 'bl3', 'avg']

# Load data
q_path = '/tf/primo/data/experiments/retrieved/enrichments'
bl_path = '/tf/primo/data/experiments/baseline/enrichments'
bl_enrichments = {}
exp_enrichments = {}

for k in queries.keys():
    exp_enrichments[k] = pd.read_csv(f'{q_path}/{k}_new.csv')
    
for k in bls:
    bl_enrichments[k] = pd.read_csv(f'{bl_path}/{k}_new.csv')


# In[4]:


# Filter out low read count sequences
seqs_to_rmv = [456, 119, 227]

bl_enrichments_rmv = {}
exp_enrichments_rmv = {}

for k,v in bl_enrichments.items():
    bl_enrichments_rmv[k] = v.loc[~v['seq_id'].isin(seqs_to_rmv)]
    
for k,v in exp_enrichments.items():
    exp_enrichments_rmv[k] = v.loc[~v['seq_id'].isin(seqs_to_rmv)]


# # Compare enrichment to predicted activity

# In[5]:


plt.figure()

for k, v in exp_enrichments_rmv.items():
    v['relative_enrichment_log'] = np.log(v['relative_enrichment'] + 1e-4)
    v['pred_activity_log'] = np.log(v['pred_activity'] + 1e-4)

    pearson_corr = v[['pred_activity_log', 'relative_enrichment_log']].corr(method='pearson').iloc[0, 1]
    spearman_corr = v[['pred_activity_log', 'relative_enrichment_log']].corr(method='spearman').iloc[0, 1]
    slope, intercept, r_value, p_value, std_err = ss.linregress(v['pred_activity_log'], v['relative_enrichment_log'])

    # Plotting
    plt.scatter(v['pred_activity_log'], v['relative_enrichment_log'], label=f'{k} (r={pearson_corr:.2f})', alpha=0.75)
    plt.plot(v['pred_activity_log'], intercept + slope * v['pred_activity_log'])
    print(spearman_corr)

plt.xlabel('Predicted Activity (Log)')
plt.ylabel('Enrichment Score (Log)')
plt.title('Predicted Activity vs. Enrichment Score')
plt.legend()

#plt.savefig('pred_v_es.svg', transparent=True , dpi=600, bbox_inches='tight', format='svg')


# # Plot heatmap

# In[6]:


df = pd.DataFrame()
num_2_seq = pd.read_csv('dfs.csv')
n2sf = num_2_seq[~num_2_seq['Unnamed: 0'].isin(seqs_to_rmv)]

df['bl1_reads'] = bl_enrichments['bl1']['reads']
df['bl2_reads'] = bl_enrichments['bl2']['reads']
for k,v in exp_enrichments_rmv.items():
    df[k+'_pred_activity'] = v['pred_activity']
    df[k+'_enr_score'] = v['relative_enrichment'].fillna(0)
    df[k+'_reads'] = v['reads']


dfm = n2sf.merge(df, left_index=True, right_index=True)


# In[7]:


yd = pd.read_csv('/tf/primo/data/simulation/targets/yields_dists_new.csv').set_index("Photo_ID")
y1 = pd.read_hdf('/tf/primo/data/targets/query_target_dists_new.h5').rename(columns = {'index': 'Photo_ID'})
yd['lego'] = y1['luis_lego']
label_lkup = {'callie_janelle' :'cat',
              'webb_deep' :'webb',
              'dna_1' :'dna',
              'frame_352' :'bigfoot',
             'yuan_taipei': 'building',
              
             'FeatureSequence': 'feature_seq'}
yd = yd.rename(columns=label_lkup)
yd['Photo_ID'] = yd.index
ydk = yd[['Photo_ID','feature_seq', 'cat', 'webb', 'dna', 'building', 'bigfoot', 'lego' ]]


# In[8]:


dfm


# In[9]:


cdf = pd.DataFrame()
cdf['feature_seq'] = dfm['feature_seq']
for k,v in exp_enrichments_rmv.items():
    cdf[k+'_enr_score'] = dfm[k+'_enr_score']
    cdf[k+'_reads'] = dfm[k+'_reads']
    cdf[k+'_pred_activity'] = dfm[k+'_pred_activity']

merged_df = pd.merge(cdf, ydk, on='feature_seq', how='left')


# In[10]:


from sklearn.preprocessing import MinMaxScaler


# Initialize the MinMaxScaler
scaler = MinMaxScaler()


# In[11]:


merged_df = merged_df.rename(columns={'bigfoot': 'frame'})


# In[17]:


merged_df


# In[18]:


queries = {'cat' : 'ACCGGTAAGGCACAGAAACG',
           'frame' : 'ATACGCAAGGAACAAAAACG'}


# In[23]:


color = ['darkblue','green']


# In[25]:


import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20, 'svg.fonttype' : 'none', 
                     'figure.figsize' : [15,8], 'xtick.major.size' : 10, 
                     'ytick.major.size' : 10})

num_bins = 200
df_sample = merged_df
for i, k in enumerate(queries.keys()):
    if k == "webb":
        continue
    #fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14, 8 * nrows))
    data = df_sample[k]
    data_filtered = data
    if not data_filtered.empty:
        plt.hist(data_filtered, bins=num_bins,label=k, alpha=0.25, color=color[i])
        q1, q2, q3 = np.percentile(data_filtered, [5, 50, 95])
        min_val = 40
        max_val = 200
    else:
        print(f'No data below 1 for {k}_enr_score')
plt.axvline(75, color='r', linestyle='--', label='Similarity Threshold')
plt.xticks(np.arange(min_val, max_val, step=20))
plt.xlim(min_val, max_val) 

plt.xticks(rotation=45)
plt.title('Image Euclidean Distances from Queries')
plt.xlabel('Euclidean Distance')
plt.ylabel('Frequency')
#plt.savefig(f"{k}_Ecu_hist.svg", transparent=True , dpi=600, bbox_inches='tight', format='svg')
plt.legend()
plt.axvline(75, color='r', linestyle='--', label='Similarity Threshold')
plt.savefig(f"Ecu_hist.svg", transparent=True , dpi=600, bbox_inches='tight', format='svg')

plt.show()


# In[30]:


'''
Heatmap that shows the percent of each bin by euclid dist that fall into each bin of enrichment score

'''
nrows = math.ceil(math.sqrt(len(exp_enrichments_rmv.keys())))
ncols = math.ceil(len(exp_enrichments_rmv.keys()) / nrows)

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15,15 * nrows))
axs = axs.ravel()

df_sample = merged_df
colors = ['Blues', 'Greens', 'Purples', 'Greens', 'Purples', 'Reds']
xlabs = ['$\leq 10^{-5}$', '$(10^{-4}, 10^{-3}]$', '$(10^{-3}, 10^{-2}]$', '$(10^{-2}, 10^{-1}]$', '$(10^{-1}, 10^{0}]$']
for i, k in enumerate(exp_enrichments_rmv.keys()):
    step_size = 2
    min_val = (df_sample[k].min() // step_size) * step_size
    #max_distance_threshold = 200
    #bins_cat = np.arange(min_val, max_distance_threshold, step=step_size)
    bins_cat = np.arange(min_val, df_sample[k].max(), step=step_size)
    #bins_cat = np.append(bins_cat, df_sample[k].max())
    df_sample[k+'_binned'] = pd.cut(df_sample[k], bins=bins_cat)
    
    bin_position_70 = np.digitize(75, bins_cat)
    print(df_sample[k].max())
    scaler = MinMaxScaler()


    df_sample[k+'_pred_activity_norm'] = scaler.fit_transform(df_sample[[k+'_pred_activity']])
    bins_enr_score = [0, .0001,.001, .01, .1, 1]
    df_sample[k+'_pred_activity_binned'] = pd.cut(df_sample[k+'_pred_activity_norm'], bins=bins_enr_score)

    hist_data = df_sample.groupby([k+'_binned', k+'_pred_activity_binned']).size().unstack(fill_value=0)
    norm_hist_data = hist_data.div(hist_data.sum(axis=1), axis=0).fillna(0)
    y_coord_70 = bin_position_70
    sns.heatmap(norm_hist_data, cmap=colors[i], ax=axs[i])
    axs[i].set_title(f'{k}')
    axs[i].set_xlabel('Predicted Activity')
    axs[i].set_ylabel('Euclidean Distance')
    
    yticks = axs[i].get_yticks()
    xticks = axs[i].get_xticks()
    current_labels = [str(item.get_text()) for item in axs[i].get_yticklabels()]
    formatted_labels = [f"({int(float(label.split(',')[0][1:]))}, {int(float(label.split(',')[1][:-1]))}]"
                        for label in current_labels]
    axs[i].set_yticklabels(formatted_labels)
    labs = axs[i].get_yticklabels()[::6]
    axs[i].set_yticks(yticks[::6])
    axs[i].set_yticklabels(labs)
    axs[i].axhline(y=y_coord_70, color='black', linestyle='--')
    axs[i].axvline(x=xticks[3]+.5, color='black', linestyle='--')
    axs[i].tick_params(axis='both', which='major', length=10)
    axs[i].tick_params(axis='both', which='major')
    axs[i].set_xticklabels(xlabs, rotation=45)
    
    for tick in xticks:
        axs[i].axvline(x=tick+.5, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

    for _, spine in axs[i].spines.items():
        spine.set_visible(True)


plt.tight_layout()
#fig.savefig('heatmap_fig6_pred.svg', transparent=True , dpi=600, bbox_inches='tight', format='svg')


# In[11]:


'''
Heatmap that shows the percent of each bin by euclid dist that fall into each bin of enrichment score

'''

nrows = math.ceil(math.sqrt(len(exp_enrichments_rmv.keys())))
ncols = math.ceil(len(exp_enrichments_rmv.keys()) / nrows)

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10 * nrows))
axs = axs.ravel()  # Flatten the axes array

df_sample = merged_df
colors = ['Blues', 'Greens', 'Purples', 'Greens', 'Purples']
xlabs = ['$\leq 10^{-5}$', '$(10^{-4}, 10^{-3}]$', '$(10^{-3}, 10^{-2}]$', '$(10^{-2}, 10^{-1}]$', '$(10^{-1}, 10^{0}]$']
for i, k in enumerate(exp_enrichments_rmv.keys()):
    step_size = 2
    min_val = (df_sample[k].min() // step_size) * step_size
    #max_distance_threshold = 200
    #bins_cat = np.arange(min_val, max_distance_threshold, step=step_size)
    bins_cat = np.arange(min_val, df_sample[k].max(), step=step_size)
    #bins_cat = np.append(bins_cat, df_sample[k].max())
    df_sample[k+'_binned'] = pd.cut(df_sample[k], bins=bins_cat)

    bin_position_70 = np.digitize(75, bins_cat)
    print(df_sample[k].max())
    scaler = MinMaxScaler()
    df_sample[k+'_enr_score_norm'] = scaler.fit_transform(df_sample[[k+'_enr_score']])
    bins_enr_score = [0,.0001, .001, .01, .1, 1]
    df_sample[k+'_enr_score_binned'] = pd.cut(df_sample[k+'_enr_score_norm'], bins=bins_enr_score)
    
    hist_data = df_sample.groupby([k+'_binned', k+'_enr_score_binned']).size().unstack(fill_value=0)
    norm_hist_data = hist_data.div(hist_data.sum(axis=1), axis=0).fillna(0)
    y_coord_70 = bin_position_70
    sns.heatmap(norm_hist_data, cmap=colors[i], ax=axs[i])
    axs[i].set_title(f'{k}')
    axs[i].set_xlabel('Enrichment Score')
    axs[i].set_ylabel('Euclidean Distance')
    xticks = axs[i].get_xticks()
    yticks = axs[i].get_yticks()
    current_labels = [str(item.get_text()) for item in axs[i].get_yticklabels()]
    formatted_labels = [f"({int(float(label.split(',')[0][1:]))}, {int(float(label.split(',')[1][:-1]))}]"
                        for label in current_labels]
    axs[i].set_yticklabels(formatted_labels)
    labs = axs[i].get_yticklabels()[::6]
    axs[i].set_yticks(yticks[::6])
    axs[i].set_yticklabels(labs)
    axs[i].axhline(y=y_coord_70, color='black', linestyle='--')
    axs[i].axvline(x=xticks[3]+.5, color='black', linestyle='--')
    axs[i].tick_params(axis='both', which='major', length=10)
    axs[i].tick_params(axis='both', which='major')
    axs[i].set_xticklabels(xlabs, rotation=45)
    for tick in xticks:
        axs[i].axvline(x=tick+.5, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)


    # Add a border around the plot
    for _, spine in axs[i].spines.items():
        spine.set_visible(True)


plt.tight_layout()
#fig.savefig('heatmap_fig6.svg', transparent=True , dpi=600, bbox_inches='tight', format='svg')


# In[12]:


# Import 3rd baseline as 'random' comparison
bl_comp = pd.read_csv(q_path + '/bl3.csv')
bl_comp = bl_comp[~bl_comp['seq_id'].isin(seqs_to_rmv)]
bl_comp['cat_pred_activity'] = dfm['cat_pred_activity']
bl_comp['bigfoot_pred_activity'] = dfm['bigfoot_pred_activity']
# Ugly, but in the sake of time good enough
bl_comp['cat_enr_score'] = bl_comp['relative_enrichment']
bl_comp['bigfoot_enr_score'] = bl_comp['relative_enrichment']


# In[13]:


fig, axs = plt.subplots(4, 1, figsize=(20, 15))


def display_colormap(cmap, data, ax, k):
    gradient = np.vstack((data[k+'_pred_activity'], data[k+'_pred_activity']))
    ax.imshow(gradient, aspect='auto', cmap=cmap, norm=matplotlib.colors.LogNorm())
    ax.set_axis_on()
    x_ticks = [0.0001, 0.001, 0.01, 0.1, .5, 1, 2]
    x_ticks.append(data[k+'_enr_score'].max())
    tick_indices = [len(data[data[k+'_enr_score'] < x]) for x in x_ticks]
    
    ax.set_xticks(tick_indices)
    ax.tick_params(axis='x', length=15, width=2)
    ax.set_xticklabels(x_ticks)
    ax.set_yticks([])

    ax.patch.set_edgecolor('black') 
    ax.patch.set_linewidth(10)
    
    spearman_corr = data[[k+'_enr_score', k+'_pred_activity']].corr(method='spearman').iloc[0, 1]
    print(spearman_corr)

tmp = dfm.sort_values(by='cat_enr_score')
ax = axs[0]
display_colormap('viridis', tmp, ax, 'cat')

ax = axs[1]
tmp = bl_comp.sort_values(by='relative_enrichment')
display_colormap('viridis', tmp, ax, 'cat')

tmp = dfm.sort_values(by='bigfoot_enr_score')
ax = axs[2]
display_colormap('viridis', tmp, ax, 'bigfoot')

ax = axs[3]
tmp = bl_comp.sort_values(by='relative_enrichment')
display_colormap('viridis', tmp, ax, 'bigfoot')


#plt.savefig(f"esvpred_bl.svg", transparent=True , dpi=600, bbox_inches='tight', format='svg')


# In[31]:


queries


# In[32]:


# Split categories for totaling
for k in queries.keys():
    merged_df[k+'_pos_or_neg'] = np.where(merged_df[k] > 75, 0, 1)
    merged_df[k+'_sim_pos_or_neg'] = np.where(merged_df[k+'_pred_activity'] <= .1, 0, 1)
    merged_df[k+'_exp_pos_or_neg'] = np.where(merged_df[k+'_enr_score'] <= 1, 0, 1)


# In[33]:


def category_totals(y_true, y_pred):
    true_positive = sum(1 for true, pred in zip(y_true, y_pred) if true == pred == 1)
    false_positive = sum(1 for true, pred in zip(y_true, y_pred) if true == 0 and pred == 1)
    false_negative = sum(1 for true, pred in zip(y_true, y_pred) if true == 1 and pred == 0)
    true_negative = sum(1 for true, pred in zip(y_true, y_pred) if true == pred == 0)
    
    precision = true_positive / (true_positive + false_positive) if true_positive + false_positive else 0
    recall = true_positive / (true_positive + false_negative) if true_positive + false_negative else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) else 0
    totalpos = sum(1 for true in y_true if true == 1)
    totalneg = sum(1 for true in y_true if true == 0)
    return [true_positive/totalpos, false_negative/totalpos, true_negative/totalneg, false_positive/totalneg]
    #return f1



# In[36]:


def category_totals(y_true, y_pred):
    true_positive = sum(1 for true, pred in zip(y_true, y_pred) if true == pred == 1)
    false_positive = sum(1 for true, pred in zip(y_true, y_pred) if true == 0 and pred == 1)
    false_negative = sum(1 for true, pred in zip(y_true, y_pred) if true == 1 and pred == 0)
    true_negative = sum(1 for true, pred in zip(y_true, y_pred) if true == pred == 0)
    
    precision = true_positive / (true_positive + false_positive) if true_positive + false_positive else 0
    recall = true_positive / (true_positive + false_negative) if true_positive + false_negative else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) else 0
    totalpos = sum(1 for true in y_true if true == 1)
    totalneg = sum(1 for true in y_true if true == 0)
    return [true_positive, false_negative, true_negative, false_positive]
    #return f1


# In[ ]:




