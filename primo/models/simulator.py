from os import read
import multiprocessing
# import pandas as pd
import numpy as np

from .cas9 import crispr_specificity

class Simulator:
    """
    Wrapper for running context (e.g. GPU, remote-execution).
    """
    def simulate(self, feature_seq_pairs):
        """
        Takes a batch of pairs of feature sequences as a pandas dataframe, and
        estimates the CAS9 specificity for each

        Input sequences should be 20nt, or a multiple of 20nt for multi-site
        evaluation.
        """
        seqlen = len(feature_seq_pairs.target_features.iloc[0])
        if seqlen > 20:
            sites = int(seqlen / 20)
            # Handle multi-site sequences
            site_scores = np.array([
                [
                    crispr_specificity(
                        p.target_features[i*20:(i+1)*20],
                        p.query_features[i*20:(i+1)*20]
                    )
                    for i in range(sites)
                ] for _, p in feature_seq_pairs.iterrows()
            ])
            # Compute combined probability
            return 1 - np.prod(1 - site_scores, axis=1)
        else:
            return np.array(
                [crispr_specificity(p.target_features, p.query_features) for _, p in feature_seq_pairs.iterrows()]
            )