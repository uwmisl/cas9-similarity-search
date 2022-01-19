import tensorflow as tf
import tensorflow_probability as tfp
import numpy as np
import primo.tools.sequences as seqtools

subpen = np.array([-1.7449405080809126, -1.275485084790358, -1.8001827224086722, -1.9323849500279549, -1.6677722398632207, -1.6537370694565101, -1.8981469677400609, -1.0814292717607923, -1.3231152511430453, -0.99840146446464273, -1.2766126030502924, -1.073338813454068, -1.5614374592181826, -1.4737507320504855, -1.298392565410591, -1.0105000195452765, -0.43349702574711524, -0.11665543376814178, -0.17370266801790191, 0.2676084623705467, 0.051835157750172757, 0.08920809165894289, 0.075459598643889569, 0.046975071077932237]).astype(np.float32)
subtrans = np.array([[ 0.        ,  1.16616601,  0.96671383,  0.94917742],       [ 0.94076049,  0.        ,  1.18426595,  0.87129983],       [ 0.58224486,  1.11064886,  0.        ,  1.04707949],       [ 0.9633753 ,  0.98895548,  1.2293125 ,  0.        ]]).astype(np.float32)
# changes matrices to account for the nucleotide string being 'ACGT', not 'ATCG' as in seqtools.bases and the rest of the PRIMO package
finkel_bases = 'ACGT'
shift = np.array([finkel_bases.index(b) for b in seqtools.bases])
subtrans = subtrans[shift,:][:, shift]
subpen = subpen[:20][::-1]

first_time_point = 12
last_time_point = 60000
first_time_point_cleavage_rate = np.log(2)/first_time_point
last_time_point_cleavage_rate = np.log(2)/last_time_point
log10_ub = np.log10(first_time_point_cleavage_rate).astype(np.float32)
log10_lb = np.log10(last_time_point_cleavage_rate).astype(np.float32)

def bandpass_hinge(x):
    return tf.maximum(log10_lb, tf.minimum(x, log10_ub))

def log10_crispr_spec(seq_pairs):

    # ensure that sequence pairs have dimension: (batch, 2 sequences, 20 nt, 4 channels)
    seq_pairs.shape.assert_is_compatible_with([None, 2, 20, 4])
    seq_len = seq_pairs.shape[2]

    # separate first and sequences from each pair
    ref = seq_pairs[:, 0, :, :]
    obs = seq_pairs[:, 1, :, :]

    # computes the outer product of one-hot vectors at each position for each pair
    # for syntax see https://numpy.org/doc/stable/reference/generated/numpy.einsum.html
    subst_ids = tf.einsum('...i,...j->...ij', ref, obs)

    # turn 4x4 match matrix into 1x16 vector
    subst_ids_flat = tf.reshape(subst_ids, [-1, seq_len, 16])

    # multiply by scores and sum
    subst_scores = tf.reduce_max(subst_ids_flat * tf.constant(subtrans.flatten()), -1)
    # subst_scores = tf.reduce_sum(subst_ids_flat * tf.constant(subtrans.flatten()), -1)

    # compute dot product of position penalty and substitution type
    scores = tf.reduce_sum(subst_scores * tf.constant(subpen.flatten()), -1)

    # adjust to range 0 to 1

    # scores = 1 + (scores - log10_ub) / (log10_ub - log10_lb)
    # return tfp.math.clip_by_value_preserve_gradient(scores, 0.0, 1.0)

    scores = tfp.math.clip_by_value_preserve_gradient(scores, log10_lb, log10_ub)
    return scores - log10_ub
    #return 10 ** (bandpass_hinge(scores + log10_ub) - log10_ub)
    #return 10 ** (tf.minimum(scores, 0) - log10_ub)
    #return 10 ** scores

def linear_crispr_spec(seq_pairs):
    return 10 ** log10_crispr_spec(seq_pairs)

def log10_norm_crispr_spec(seq_pairs):
    return 1.0 + log10_crispr_spec(seq_pairs) / (log10_ub - log10_lb)

def crispr_spec_for_loss(seq_pairs):

    # ensure that sequence pairs have dimension: (batch, 2 sequences, 20 nt, 4 channels)
    seq_pairs.shape.assert_is_compatible_with([None, 2, 20, 4])
    seq_len = seq_pairs.shape[2]

    # separate first and sequences from each pair
    ref = seq_pairs[:, 0, :, :]
    obs = seq_pairs[:, 1, :, :]

    # computes the outer product of one-hot vectors at each position for each pair
    # for syntax see https://numpy.org/doc/stable/reference/generated/numpy.einsum.html
    subst_ids = tf.einsum('...i,...j->...ij', ref, obs)

    # turn 4x4 match matrix into 1x16 vector
    subst_ids_flat = tf.reshape(subst_ids, [-1, seq_len, 16])

    # multiply by scores and sum
    subst_scores = tf.reduce_sum(subst_ids_flat * tf.constant(subtrans.flatten()), -1)

    # compute dot product of position penalty and substitution type
    scores = tf.reduce_sum(subst_scores * tf.constant(subpen.flatten()), -1)

    return scores - log10_ub
