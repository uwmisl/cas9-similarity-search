import pandas as pd
import numpy as np
import math
from primo.models.simulator import Simulator

#Predicted Cas9 Activity
def get_activation_score(seq1, seq2):
    """
    Given two DNA sequences of length 20,
    Returns a float between 0 and 1 of the Cas9 activity between them if one were gRNA and one were the DNA seq
       Note that order doesn't matter, and 0 means minimum cleavage, 1 means max cleavage
    """
    simulator = Simulator()
    df = pd.DataFrame({
        "target_features": [seq1],
        "query_features": [seq2]
    })
    return simulator.simulate(df)[0] # this is a float 

def calc_pred_enrichment(query):
    """
    Calls get_activation_score on a given query sequence and the entire oligo pool to
    return the expected cut activity
    
    Takes a tuple or array of a query sequence and name
    """
    seq_df = pd.read_csv('/tf/primo/notebooks/03_simulation/oligos.csv')
    seq_df = seq_df.rename(columns = {'Unnamed: 0': 'Photo_ID'})
    # get all the cas sites from the seq_df dataframe
    cas_sites = seq_df['FeatureSequence'].unique()
    # make a list of cas9 activity scores
    cas9_activity = [get_activation_score(seq, query) for seq in cas_sites]

    # make list of query info to add to dataframe
    query_seq = [query]*len(cas9_activity)

    # make the dataframe with all information
    zipped = list(zip(query_seq, cas_sites, cas9_activity))
    predicted_cas9_activity_df = pd.DataFrame(zipped, columns = 
                                              ['Query_Seq','Target_Seq','wtCas9_Predicted_Activity'])
    return predicted_cas9_activity_df['wtCas9_Predicted_Activity']


# Enrichment
def calc_pool_enrichment(aligned_reads):
    """
    For each unique sequence (1-457), calculates the reads that aligned to it as a percentage of
    the total reads within the pool
    """
    aligned_reads_filter = aligned_reads[aligned_reads['alignment_genome'] != '*']
    counts = pd.to_numeric(aligned_reads_filter['alignment_genome']).value_counts(sort=False)
    result = [counts[i] if i in counts.index else 0 for i in range(457)]
    
    total_reads = sum(counts)
    pool_enrichment = counts/total_reads
    
    return counts, pool_enrichment

def calc_relative_enrichment(exp_enrichment, baseline_enrichment):
    """
    For each experiment, calculates the enrichment relative to the baseline enrichment
    """
    return exp_enrichment/baseline_enrichment.sort_index().fillna(0)

def summarize_enrichment(exp, pred_activity=None, bl=None):
    reads = calc_pool_enrichment(exp)
    df = pd.DataFrame()
    df['seq_id'] = range(457)
    df['reads'] = reads[0]
    df['pool_enrichment'] = reads[1]
    
    # if a baseline enrichment is passed in, calculate the relative enrichment
    if bl is not None:
        df['relative_enrichment'] = calc_relative_enrichment(df['pool_enrichment'], bl['pool_enrichment'])
   
    return df.fillna(0)

# Read Totals
def read_summary(alignments):
    """
    Given a guppy summary of total reads in a sequencing run, this returns the number in each of six categories:
    - Uncut +/-
    - Cut, First half +/-
    - Cut, Second half +/-
    - Unaligned
    """
    counts_df = pd.DataFrame(columns=['Experiment ID',
                                      'First, +',
                                      'First, -',
                                      'Second, +',
                                      'Second, -',
                                      'Uncut',
                                      'Unaligned',
                                      'Total Reads'])
    for k,v in alignments.items():        
        # clean up unaligned reads
        all_aligned = v[v['alignment_genome'] != '*']
        unaligned = v[v['alignment_genome'] == '*']
        
        row = [k]
        # calculates counts of each group
        fp = all_aligned[(all_aligned['alignment_genome_end'] < 56) & (all_aligned['alignment_direction'] == '+')]
        fm = all_aligned[(all_aligned['alignment_genome_end'] < 56) & (all_aligned['alignment_direction'] == '-')]
        sp = all_aligned[(all_aligned['alignment_genome_start'] >= 56) & (all_aligned['alignment_direction'] == '+')]
        sm = all_aligned[(all_aligned['alignment_genome_start'] >= 56) & (all_aligned['alignment_direction'] == '-')]
        up = all_aligned[(all_aligned['alignment_genome_start'] < 56) & (all_aligned['alignment_genome_end'] >= 56)]
       

        row.append(len(fp))
        row.append(len(fm))
        row.append(len(sp))
        row.append(len(sm))
        row.append(len(up))
        row.append(len(unaligned))
        row.append(len(all_aligned) + len(unaligned))
        # Add to data frame
        counts_df.loc[len(counts_df)] = row
    return counts_df

# calculate distance from y=x
def dist(point, a, b, c):
    return (a * point[0] + b * point[1] + c) / (math.sqrt(a * a + b * b))

# Calculate Kendall tau distance
def normalised_kendall_tau_distance(values1, values2):
    """Compute the Kendall tau distance."""
    n = len(values1)
    assert len(values2) == n, "Both lists have to be of equal length"
    i, j = np.meshgrid(np.arange(n), np.arange(n))
    a = np.argsort(values1)
    b = np.argsort(values2)
    ndisordered = np.logical_or(np.logical_and(a[i] < a[j], b[i] > b[j]), np.logical_and(a[i] > a[j], b[i] < b[j])).sum()
    return ndisordered / (n * (n - 1))


def kt_dist(plot_df):
    plot_df['x_dist'] = [dist(x, -1, 1, 0) for x in zip(plot_df.x_vals, plot_df.y_vals)]
    plot_df['enrichment_score'] = plot_df['y_vals']/plot_df['x_vals']
 
    #Sort by distance where larger distance is first and extract seq id order
    dist_rank = plot_df.sort_values(by='x_dist', ascending=False)['seq_id'].to_numpy()
    # Get rank when ordering by pred cas9 activity
    pred_rank = plot_df.sort_values(by='color',ascending=False)['seq_id'].to_numpy()
    # Get rank by enrichment score
    es_rank = plot_df.sort_values(by='enrichment_score',ascending=False)['seq_id'].to_numpy()
    return normalised_kendall_tau_distance(dist_rank, pred_rank)