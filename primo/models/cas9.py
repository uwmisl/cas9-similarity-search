import pickle
import numpy as np
from pathlib import Path

def model_values_path():
    """Returns the location of the cas9 cleavage model values.

    Returns
    -------
    os.pathlike
        Path to the pickle containing the cleavage values.
    """
    THIS_DIR = Path(__file__).parent
    path_to_model_values = Path(THIS_DIR, "cas9_cleavage_model_values.pickle")
    return path_to_model_values

# The functions contained here come from the paper: https://www.nature.com/articles/s41587-020-0646-5
# Specifically from /nucleaseq/modeling.py

# This file allows you to pass in the sequences to be compared to each other and
# outputs a cleavage score.

# Length (in seconds) of the shortest cas9 cleavage experiment.
first_time_point = 12
# Length (in seconds) of the longest cas9 cleavage experiment.
last_time_point = 60_000
first_time_point_cleavage_rate = np.log(2)/first_time_point
last_time_point_cleavage_rate = np.log(2)/last_time_point
log10_ub = np.log10(first_time_point_cleavage_rate)
log10_lb = np.log10(last_time_point_cleavage_rate)


def build_cro(ref_seq, obs_seq):
    """
    This function comes from the paper: https://www.nature.com/articles/s41587-020-0646-5
    """
    bases_and_deletion = 'ACGT-'
    bases = 'ACGT'
    sub_cro, del_cro, ins_cro = [], [], []
    i = 0
    for rb, ob in zip(ref_seq, obs_seq):
        if rb != ob:
            ri = bases_and_deletion.index(rb)
            oi = bases_and_deletion.index(ob)
            cro = (i, ri, oi)
            if rb in bases and ob in bases:
                sub_cro.append(cro)
            elif rb == '-':
                ins_cro.append(cro)
            elif ob == '-':
                del_cro.append(cro)
            else:
                raise ValueError('Unexpected input: ({}, {}, {})'.format(i, rb, ob))
        if rb != '-':
            i += 1
    return sub_cro, del_cro, ins_cro

def bandpass_hinge(x):
    return max(log10_lb, min(x, log10_ub))

def single_effects(pam_cro, sub_cro, del_cro, ins_cro,
                pampwm, subpen, subtrans, delpen, inspen, insweight):
    """
    This function comes from the paper: https://www.nature.com/articles/s41587-020-0646-5
    """
    score = log10_ub
    for i, ri, oi in pam_cro:
        score += pampwm[oi][i]
    for i, ri, oi in sub_cro:
        score += subpen[i] * subtrans[ri, oi]
    for i, ri, oi in del_cro:
        score += delpen[i]
    for i, ri, oi in ins_cro:
        score += inspen[i] * insweight[oi]
    return bandpass_hinge(score) - log10_ub

def crispr_specificity( grna_seq, dna_seq, cas_protein='wt', pam_sequence="tgg"):
    """
    This function comes from the paper: https://www.nature.com/articles/s41587-020-0646-5
    The function was originally named "log10_crisp_specificity" in their work.


    The main model function.

    Input:
        cas_protein      :str: The protein to model. For now, only wild-type ('wt') is supported, original paper includes: WT, Enh, Hypa, HF1, Cas12a
        pam_sequence     :str: The PAM sequence with not indels.
        grna_seq    :str: The aligned, 5'->3' guide RNA sequence, with insertions as hyphens
        dna_seq :str: The aligned, 5'->3' DNA target sequence (no PAM), with deletions as hyphens
    """

    # input cleaning

    # ensures that PAM sequence is uppercase and gRNA sequence is uppercase and
    # if a U is given it's now a T
    pam_sequence = pam_sequence.upper()
    grna_seq = grna_seq.upper().replace('U', 'T')
    dna_seq = dna_seq.upper()

    # specifies the protein in use (adapted from https://www.nature.com/articles/s41587-020-0646-5)
    # currently only wt cas9 is supported
    cas_protein = cas_protein.lower()
    if cas_protein != 'wt':
        raise ValueError('This protein is not currently supported. See simulator.py for instructions.')
        # if you want to incorporate more proteins, see https://www.nature.com/articles/s41587-020-0646-5
    if cas_protein == 'wt':
        # reverse the sequences (NOT compliment)
        grna_seq = grna_seq[::-1]
        dna_seq = dna_seq[::-1]

    # loads in the wt Cas9 cleavage model data
    cleavage_data_path = model_values_path()
    with open(cleavage_data_path, 'rb') as f:
        ref_param = pickle.load(f, encoding="bytes")
    params = ref_param

    # checks the PAM sequence is valid
    if pam_sequence != 'TGG':
        raise ValueError('This PAM sequence is not supported. We currently only support "TGG".')
    ref_pam = pam_sequence

    # further checks PAM
    pam_cro, pam_del_cro, pam_ins_cro = build_cro(ref_pam, pam_sequence)
    if pam_del_cro or pam_ins_cro:
        raise ValueError('Invalid PAM: {}'.format(pam_sequence))

    # Starts analyzing the sequences to be assigned a score.
    sub_cro, del_cro, ins_cro = build_cro(grna_seq, dna_seq)
    log_score = single_effects(pam_cro, sub_cro, del_cro, ins_cro, *params)
    score = 10**log_score
    return score


if __name__ == "__main__":
    grna_seq = "ACCGCCCCATCGACCTAGCA"
    dna_seq = "ATCGCCCCATCGACCTAGCA"
    score = crispr_specificity(grna_seq, dna_seq)
    print(f"score: {score}")