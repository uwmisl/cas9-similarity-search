#!/usr/bin/env python
# coding: utf-8

# In[18]:


import pandas as pd
import numpy as np
from datetime import date
import os
from tqdm import tqdm # for progress bar
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import random 
from Bio import pairwise2
import csv
import sys
from sklearn.linear_model import LinearRegression as LR
import statsmodels.api as sm
import scipy.stats as ss
import pylab
import mplcursors
import math
from IPython.core.display import display, HTML
import io
import zipfile
from PIL import Image
import pickle
from data.sequencing import data_analysis
from data.sequencing import visualization
from data.sequencing import data_import
import matplotlib.patches as patches


# # Import non-seq data

# In[19]:


dists = data_import.import_dists()
dfs = pd.read_csv('dfs.csv')
dfs


# In[34]:


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


# In[21]:


# Filter out low read count sequences
seqs_to_rmv = [456, 119, 227]

bl_enrichments_rmv = {}
exp_enrichments_rmv = {}

for k,v in bl_enrichments.items():
    bl_enrichments_rmv[k] = v.loc[~v['seq_id'].isin(seqs_to_rmv)]
    
for k,v in exp_enrichments.items():
    exp_enrichments_rmv[k] = v.loc[~v['seq_id'].isin(seqs_to_rmv)]


# # Compare enrichment to predicted activity

# In[22]:


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


# # Plot heatmap

# In[23]:


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


# In[24]:


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


# In[25]:


new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none',
}
merged_df = merged_df.rename(columns={'bigfoot': 'frame'})


# In[26]:


cdf = pd.DataFrame()
cdf['feature_seq'] = dfm['feature_seq']
for k,v in exp_enrichments_rmv.items():
    cdf[k+'_enr_score'] = dfm[k+'_enr_score']
    cdf[k+'_reads'] = dfm[k+'_reads']
    cdf[k+'_pred_activity'] = dfm[k+'_pred_activity']

merged_df = pd.merge(cdf, ydk, on='feature_seq', how='left')


# In[29]:


merged_df = merged_df.rename(columns={'bigfoot': 'frame'})


# In[30]:


for k in queries.keys():
    merged_df[k+'_pos_or_neg'] = np.where(merged_df[k] > 75, 0, 1)
    merged_df[k+'_sim_pos_or_neg'] = np.where(merged_df[k+'_pred_activity'] < .2, 0, 1)
    merged_df[k+'_exp_pos_or_neg'] = np.where(merged_df[k+'_enr_score'] <= 1, 0, 1)
    
    
    


# In[31]:


def manual_f1_score(y_true, y_pred):
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



# In[35]:


for k in queries.keys():
    score_manual_exp = manual_f1_score(merged_df[k+'_pos_or_neg'], merged_df[k+'_exp_pos_or_neg'])
    
    print(k+':')
    print('%.6f' % score_manual_exp[0])
    print('%.6f' % score_manual_exp[1])
    print('%.6f' % score_manual_exp[2])
    print('%.6f' % score_manual_exp[3])
    print()
    print()


# In[36]:


for k in queries.keys():
    score_manual_exp = manual_f1_score(merged_df[k+'_pos_or_neg'], merged_df[k+'_exp_pos_or_neg'])
    
    print(k+':')
    print('%.6f' % score_manual_exp[0])
    print('%.6f' % score_manual_exp[1])
    print('%.6f' % score_manual_exp[2])
    print('%.6f' % score_manual_exp[3])
    print()
    print()


# In[38]:


merged_df


# In[39]:


merged_df.value_counts('feature_seq').index


# In[40]:


df = pd.DataFrame()
df['feature_seq'] = merged_df.value_counts('feature_seq').index
df['count'] = merged_df.value_counts('feature_seq').values

for k in queries.keys():
    tmp_sim = []
    tmp_diff = []
    min_l = []
    max_l = []
    es_l = []
    pa_l = []
    print(k)
    for i in merged_df.value_counts('feature_seq').index:
        tmp = merged_df[merged_df['feature_seq'] == i]
        num_sim = len(tmp[tmp[k] <= 75])
        num_diff = len(tmp[tmp[k] > 75])
        min_ = tmp[k].min()
        max_ = tmp[k].max()
        es = tmp[k+'_enr_score'].min()
        pa = tmp[k+'_pred_activity'].max()
        
        es_l.append(es)
        pa_l.append(pa)
        min_l.append(min_)
        max_l.append(max_)
        tmp_sim.append(num_sim)
        tmp_diff.append(num_diff)
    
    df[k+'_pa'] = pa_l
    df[k+'_es'] = es_l
    df[k+'_min'] = min_l
    df[k+'_max'] = max_l
    df[k+'_sim'] = tmp_sim
    df[k+'_diff'] = tmp_diff
    
        


# In[41]:


dfm


# In[42]:


newdf = pd.merge(dfm, df, on = 'feature_seq')


# In[43]:


newdf


# In[ ]:


newdf.to_csv('bin_summary_2.csv')


# In[ ]:


df_sample


# In[ ]:


merged_df


# In[ ]:


for k in queries.keys():
    recall_data = []
    false_positive_data = []
    thresholds = merged_df[k+'_pred_activity'].unique()
    thresholds.sort()
    pairs_s = merged_df[merged_df[k] <= 75]
    pairs_d = merged_df[merged_df[k] > 75]

    for t in thresholds:
        recall = len(pairs_s[pairs_s[k+'_pred_activity'] > t])/len(pairs_s)
        false_positive = len(pairs_d[pairs_d[k+'_pred_activity'] > t])/len(pairs_d)
        recall_data.append(recall)
        false_positive_data.append(false_positive)
    plt.figure(figsize=(10, 5))   
    plt.scatter(false_positive_data, recall_data, c=thresholds, cmap=plt.cm.viridis_r)
    plt.plot(false_positive_data, false_positive_data, linestyle='dashed', color='grey')
    plt.ylabel("Proportion Correctly Retrieved")
    plt.xlabel("Proportion Mistakenly Retrieved")
    plt.colorbar(label="Yield Threshold")
    plt.title("20nt Cas9 Model - "+k)


# In[ ]:


merged_df


# In[44]:


for k in queries.keys():
    merged_df[k] = merged_df[k].astype(int)
    print(len(merged_df))
    newdf = merged_df.drop_duplicates(subset=[k, k+'_pred_activity'])
    print(len(newdf))
    enrs = newdf[(newdf[k+'_pred_activity'] > .1) & (newdf[k] <= 75)].sort_values(by=k)
    deps = newdf[(newdf[k+'_pred_activity'] <= .1) & (newdf[k] <= 75)].sort_values(by=k)
    enrd = newdf[(newdf[k+'_pred_activity'] > .1) & (newdf[k] > 75)].sort_values(by=k)
    depd = newdf[(newdf[k+'_pred_activity'] <= .1) & (newdf[k] > 75)].sort_values(by=k)

    plt.figure(figsize=(10, 10)) 
    plt.scatter(enrs[k], enrs[k+'_pred_activity'], label=k+' enrs', s=20)
    plt.scatter(deps[k], deps[k+'_pred_activity'], label=k+' deps', s=20)
    plt.scatter(enrd[k], enrd[k+'_pred_activity'], label=k+' enrd', s=20)
    plt.scatter(depd[k], depd[k+'_pred_activity'], label=k+' depd', s=20)

    plt.axvline(x=75, color='b', linestyle='-')
    plt.axhline(y=.1, color='b', linestyle='-')
    plt.legend()
    plt.rcParams.update(new_rc_params)
    #plt.savefig(k+".svg", transparent=True, dpi=600, bbox_inches='tight', format='svg')
    # Clear the set for the next iteration
    


# In[45]:


merged_df['cat_pred_activity'].value_counts()


# In[46]:


value_counts.index


# In[47]:


import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

for k in queries.keys():
    merged_df[k] = merged_df[k].astype(int)
    newdf = merged_df.drop_duplicates(subset=[k, k+'_pred_activity'])
    enrs = newdf[(newdf[k+'_pred_activity'] > .1) & (newdf[k] <= 75)].sort_values(by=k)
    deps = newdf[(newdf[k+'_pred_activity'] <= .1) & (newdf[k] <= 75)].sort_values(by=k)
    enrd = newdf[(newdf[k+'_pred_activity'] > .1) & (newdf[k] > 75)].sort_values(by=k)
    depd = newdf[(newdf[k+'_pred_activity'] <= .1) & (newdf[k] > 75)].sort_values(by=k)

    fig, (ax_scatter, ax_hist) = plt.subplots(1, 2, gridspec_kw={"width_ratios": [5, 1]}, figsize=(10, 8))

    ax_scatter.scatter(enrs[k], enrs[k+'_pred_activity'], label=k+' enrs', s=20)
    ax_scatter.scatter(deps[k], deps[k+'_pred_activity'], label=k+' deps', s=20)
    ax_scatter.scatter(enrd[k], enrd[k+'_pred_activity'], label=k+' enrd', s=20)
    ax_scatter.scatter(depd[k], depd[k+'_pred_activity'], label=k+' depd', s=20)
    ax_scatter.axvline(x=75, color='black', linestyle='--')
    ax_scatter.axhline(y=.1, color='black', linestyle='--')
    ax_scatter.legend()

    min_x_limit = 0.1

    ax_hist.hist(merged_df[k+'_pred_activity'], bins=40, orientation='horizontal')
    ax_hist.set_xscale('log')
    ax_hist.set_ylim(ax_scatter.get_ylim())
    ax_hist.set_xlim(left=min_x_limit)  # Set the lower limit for the x-axis
    ax_hist.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=5))
    ax_hist.grid(which='major', linestyle='--')

    ##plt.savefig(k+".svg", transparent=True, dpi=600, bbox_inches='tight', format='svg')
    plt.show()


# In[ ]:


for k in queries.keys():
    num_enr = len(merged_df[merged_df[k+'_enr_score'] >1])
    num_sim = len(merged_df[merged_df[k] <75])
    print(k)
    print(num_enr)
    print(num_sim)


# In[ ]:


df_sample


# In[ ]:


for k in queries.keys():
    print(len(df_sample[(df_sample[k+'_pred_activity_norm'] > .1) & (df_sample[k] <= 75)]))
    print(len(df_sample[(df_sample[k+'_pred_activity_norm'] <= .1) & (df_sample[k] <= 75)]))
    print(len(df_sample[(df_sample[k+'_pred_activity_norm'] > .1) & (df_sample[k] > 75)]))
    print(len(df_sample[(df_sample[k+'_pred_activity_norm'] <= .1) & (df_sample[k] > 75)]))
    
    print()
    
    print(len(df_sample[(df_sample[k+'_enr_score_norm'] > .1) & (df_sample[k] <= 75)]))
    print(len(df_sample[(df_sample[k+'_enr_score_norm'] <= .1) & (df_sample[k] <= 75)]))
    print(len(df_sample[(df_sample[k+'_enr_score_norm'] > .1) & (df_sample[k] > 75)]))
    print(len(df_sample[(df_sample[k+'_enr_score_norm'] <= .1) & (df_sample[k] > 75)]))
    
    print()


# In[ ]:


for k in queries.keys():
    print(len(df_sample[(df_sample[k+'_pred_activity'] > .1) & (df_sample[k] <= 75)]))
    print(len(df_sample[(df_sample[k+'_pred_activity'] <= .1) & (df_sample[k] <= 75)]))
    print(len(df_sample[(df_sample[k+'_pred_activity'] > .1) & (df_sample[k] > 75)]))
    print(len(df_sample[(df_sample[k+'_pred_activity'] <= .1) & (df_sample[k] > 75)]))
    
    print()
    
    print(len(df_sample[(df_sample[k+'_enr_score'] > 1) & (df_sample[k] <= 75)]))
    print(len(df_sample[(df_sample[k+'_enr_score'] <= 1) & (df_sample[k] <= 75)]))
    print(len(df_sample[(df_sample[k+'_enr_score'] > 1) & (df_sample[k] > 75)]))
    print(len(df_sample[(df_sample[k+'_enr_score'] <= 1) & (df_sample[k] > 75)]))
    
    print()


# In[ ]:




