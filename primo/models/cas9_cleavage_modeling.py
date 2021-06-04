import pickle

# The functions contained here come from the paper: https://www.nature.com/articles/s41587-020-0646-5
# Specifically from /nucleaseq/modeling.py

# This file allows you to pass in the sequences to be compared to each other and
# outputs a cleavage score.

first_time_point = 12
last_time_point = 60000
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

def log10_crispr_specificity(prot='wt', pam_seq, grna_seq, nts_dna_seq):
    """
    This function comes from the paper: https://www.nature.com/articles/s41587-020-0646-5
    The function was originally named "log10_crisp_specificity" in their work.

    The main model function.

    Input:
        prot        :str: The protein to model. Options: WT, Enh, Hypa, HF1, Cas12a
        pam_seq     :str: The NTS PAM sequence with not indels
        grna_seq    :str: The aligned, 5'->3' guide RNA sequence, with insertions as hyphens
        nts_dna_seq :str: The aligned, 5'->3' NTS DNA target sequence (no PAM), with deletions as hyphens
    """
    # input cleaning

    # ensures that PAM sequence is uppercase and gRNA sequence is uppercase and
    # if a U is given it's now a T
    pam_seq = pam_seq.upper()
    grna_seq = grna_seq.upper().replace('U', 'T')
    nts_dna_seq = nts_dna_seq.upper()

    # specifies the protein in use (adapted from https://www.nature.com/articles/s41587-020-0646-5)
    # currently only wt cas9 is supported
    prot = prot.lower()
    if prot != 'wt':
        raise ValueError('This protein is not currently supported. See simulator.py for instructions.')
        # if you want to incorporate more proteins, see https://www.nature.com/articles/s41587-020-0646-5
    if prot == 'wt':
        # reverse the sequences (NOT compliment)
        grna_seq = grna_seq[::-1]
        nts_dna_seq = nts_dna_seq[::-1]

    # loads in the wt Cas9 cleavage model data
    with open('cas9_cleavage_model_values.pickle', 'rb') as f:
        ref_param = pickle.load(f)
    params = ref_param

    # checks the PAM sequence is valid
    if pam_seq != 'TGG':
        raise ValueError('This PAM sequence is not supported. We currently only support "TGG".')
    ref_pam = pam_seq

    # further checks PAM
    pam_cro, pam_del_cro, pam_ins_cro = build_cro(ref_pam, pam_seq)
    if pam_del_cro or pam_ins_cro:
        raise ValueError('Invalid PAM: {}'.format(pam_seq))

    # prepares the sequences to be compared and scored
    simple_rna = grna_seq.replace('-', '')
    simple_dna = nts_dna_seq.replace('-', '')
    if (not set(simple_rna + simple_dna) <= ('A', 'T', 'C','G')):
        raise ValueError('Invalid input sequences:\n{}\n{}'.format(grna_seq, nts_dna_seq))
    sub_cro, del_cro, ins_cro = build_cro(grna_seq, nts_dna_seq)
    return single_effects(pam_cro, sub_cro, del_cro, ins_cro, *params)