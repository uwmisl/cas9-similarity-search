import tensorflow as tf
import tensorflow_probability as tfp
import numpy as np
import primo.tools.sequences as seqtools

subpen = np.array([-1.7449405080809126, -1.275485084790358, -1.8001827224086722, -1.9323849500279549, -1.6677722398632207, -1.6537370694565101, -1.8981469677400609, -1.0814292717607923, -1.3231152511430453, -0.99840146446464273, -1.2766126030502924, -1.073338813454068, -1.5614374592181826, -1.4737507320504855, -1.298392565410591, -1.0105000195452765, -0.43349702574711524, -0.11665543376814178, -0.17370266801790191, 0.2676084623705467, 0.051835157750172757, 0.08920809165894289, 0.075459598643889569, 0.046975071077932237]).astype(np.float32)
subtrans = np.array([[ 0.        ,  1.16616601,  0.96671383,  0.94917742],       [ 0.94076049,  0.        ,  1.18426595,  0.87129983],       [ 0.58224486,  1.11064886,  0.        ,  1.04707949],       [ 0.9633753 ,  0.98895548,  1.2293125 ,  0.        ]])
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
data_span = log10_ub - log10_lb

def bandpass_hinge(x):
    return tf.maximum(log10_lb, tf.minimum(x, log10_ub))

def log10_crispr_spec(seq_pairs):

    # ensure that sequence pairs have dimension: (batch, 2 sequences, 20 nt, 4 channels)
    seq_pairs.shape.assert_is_compatible_with([None, 2, 20, 4]) #NEVER change shape to 80 or any other multiple of 20!
    seq_len = seq_pairs.shape[2]

    # separate first and sequences from each pair
    ref = seq_pairs[:, 0, :, :]
    obs = seq_pairs[:, 1, :, :]

    # computes the outer product of one-hot vectors at each position for each pair
    # for syntax see https://numpy.org/doc/stable/reference/generated/numpy.einsum.html
    subst_ids = tf.einsum('...i,...j->...ij', ref, obs)

    # turn 4x4 match matrix into 1x16 vector
    subst_ids_flat = tf.reshape(subst_ids, [-1, seq_len, 16])

    # multiply by scores and take highest value
    subst_scores = tf.reduce_max(subst_ids_flat * subtrans.flatten(), -1)

    # compute dot product of position penalty and substitution type
    scores = tf.reduce_sum(subst_scores * subpen.flatten(), -1)

    # threshold final result
    return bandpass_hinge(scores + log10_ub) - log10_ub

def linear_crispr_spec(mid_point=None):
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
        return power ** log10_crispr_spec(seq_pairs)

    return f

def log10_norm_crispr_spec(seq_pairs):
    logscores = log10_crispr_spec(seq_pairs)
    return tfp.math.clip_by_value_preserve_gradient(1.0 + logscores / (log10_ub - log10_lb), 0.0, 1.0)

def dotproduct_crispr_spec(seq_pairs):
    """Alternative predictor function using dot-product for base comparison

    Normalize the softmax inputs, and take the penalty at each position to be
    one minus the dot product of the two "base vectors", weighted according to the
    position dependent substitution penalty from the finkelstein paper.

    This `log10_crispr_spec` has the same output value for one-hot inputs, but
    varies for uncertain softmax inputs.

    With `log10_crispr_spec`, an input such as [0.8, 0.066, 0.066, 0.0.66] / [0.8, 0.066, 0.066, 0.0.66]
    at a single base results in a substantial penalty to cleave rate. In this
    implementation, these would be treated as matching, and not penalize the rate.
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
    scores = tf.minimum(0.0, scores)
    return scores

def dotproduct_linearized(mid_point=None):
    """Adjusted dotproduct_crispr_spec output

    Converts the output of `dotproduct_crispr_spec` to be on range 0 to 1.

    This weights the output by confidence. At low confidence the output will be
    0.5 -- ambiguous. As the confidence of the input bases increases, it
    adjusts the output towards the score value. This creates a smoother loss
    function for better training.
    """
    power = 10.0
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

def make_multisite_predictor(predictor):
    """Return a function that will predict multiple sites using the predictor
    function provided to score each

    Combined probability is 1 - (1 - P(site1)) * (1 - P(site2)) ... * (1 - P(siteN))
    """
    def multisite_predict(seq_pairs):
        n_sites = int(seq_pairs.shape[2] / 20)
        # Split into separate sites
        sites = tf.stack(tf.split(seq_pairs, n_sites, axis=2))
        # Apply predictor to sites independently
        scores = tf.map_fn(tf.function(predictor), sites)
        # Compute combined probability
        return 1 - tf.reduce_prod(1 - scores, axis=0)

    return multisite_predict

def log_multisite_predictor(seq_pairs):
    n_sites = int(seq_pairs.shape[2] / 20)
    # Split into separate sites
    sites = tf.stack(tf.split(seq_pairs, n_sites, axis=2))
    # Apply predictor to sites independently
    # TODO: Try using log10_crispr_spec instead of dotproduct_crispr_spec. It
    # should be the same with hardmax inputs
    scores = tf.map_fn(tf.function(dotproduct_crispr_spec), sites)
    linear_scores = 10**scores
    return tf.experimental.numpy.log10(tf.reduce_sum(linear_scores, axis=0))
