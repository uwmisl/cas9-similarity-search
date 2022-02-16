import numpy as np
import tensorflow as tf
from ..tools import sequences as seqtools

subpen = np.array([-1.7449405080809126, -1.275485084790358, -1.8001827224086722, -1.9323849500279549, -1.6677722398632207, -1.6537370694565101, -1.8981469677400609, -1.0814292717607923, -1.3231152511430453, -0.99840146446464273, -1.2766126030502924, -1.073338813454068, -1.5614374592181826, -1.4737507320504855, -1.298392565410591, -1.0105000195452765, -0.43349702574711524, -0.11665543376814178, -0.17370266801790191, 0.2676084623705467, 0.051835157750172757, 0.08920809165894289, 0.075459598643889569, 0.046975071077932237]).astype(np.float32)
subpen = subpen[:20][::-1]
subtrans = np.array([[ 0.        ,  1.16616601,  0.96671383,  0.94917742],       [ 0.94076049,  0.        ,  1.18426595,  0.87129983],       [ 0.58224486,  1.11064886,  0.        ,  1.04707949],       [ 0.9633753 ,  0.98895548,  1.2293125 ,  0.        ]]).astype(np.float32)
# changes matrices to account for the nucleotide string being 'ACGT', not 'ATCG' as in seqtools.bases and the rest of the PRIMO package
finkel_bases = 'ACGT'
shift = np.array([finkel_bases.index(b) for b in seqtools.bases])
subtrans = subtrans[shift,:][:, shift]

def dotproduct_crispr_spec(seq_pairs):
    """Alternative cleave rate model

    Normalize the softmax inputs, and take the penalty at each position to be
    one minus the dot product of the two "base vectors", weighted according to the
    position dependent substitution penalty from the finkelstein paper.

    This ignores the base dependent weighting -- i.e. difference between substituting
    A, T, C, or G -- but allows base outputs which are less confident but still
    similar to be evaluated as such during training.

    With the previous function, an input like [0.8, 0.066, 0.066, 0.0.66] / [0.8, 0.066, 0.066, 0.0.66]
    at a single base results in a substantial penalty to cleave rate.
    """
    ref = seq_pairs[:, 0, :, :]
    obs = seq_pairs[:, 1, :, :]
    ref_norm, _ = tf.linalg.normalize(ref, axis=-1)
    obs_norm, _ = tf.linalg.normalize(obs, axis=-1)

    m1 = tf.einsum('bij,jk->bijk', ref, subtrans)
    m2 = tf.einsum('bij,bikj->bijk', obs, m1)
    base_factors = tf.reduce_sum(tf.reduce_sum(m2, axis=-1), axis=-1)
    # 1 - Dot product
    x = 1 -  tf.reduce_sum(ref_norm*obs_norm, axis=-1)

    scores = tf.reduce_sum(x * base_factors * tf.constant(subpen.flatten(), dtype=tf.float32), -1)
    scores = tf.minimum(0, scores)
    return scores

def dotproduct_linearized(mid_point=None):
    """Returns a predictor function which will scale the log10 scores such that the
    given `mid_point` value is 0.5.

    mid_point is the output of log10_crispr_spec, and shoudl be on range [log10_lb, log10_ub]

    If mid_point is none, cleave rate is linearized; i.e. return 10**log_10_crispr_spec(x)
    """
    power = 10
    if mid_point is not None:
        power = 0.5 ** (1 / mid_point)

    def f(seq_pairs):
        """
        seq_pairs batch_size x 2 x SEQLEN x 4
        """
        scores = power ** dotproduct_crispr_spec(seq_pairs)
        confidence = tf.reduce_mean(tf.reduce_max(seq_pairs, -1), axis=2)
        confidence = tf.reduce_prod(confidence, axis=1)
        return 0.5 + (scores - 0.5) * confidence

    return f