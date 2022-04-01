from os import read
import multiprocessing
import numpy as np
import pandas as pd

from .cas9 import crispr_specificity

from ..tools import sequences as seqtools

def f(row):
    return crispr_specificity(row[1].target_features, row[1].query_features)

class Simulator:
    """
    Wrapper for running context (e.g. GPU, remote-execution).
    """

    # defaults = {
    #     # Reverse Primer.
    #     "RP": "GTCCTCAACAACCTCCTG",
    #     # First 6 nucleotides of the reverse primer.
    #     "toehold": "GTCCTC",

    #     # Target molar concentration.
    #     "t_conc": 1e-9,
    #     # Query molar concentration.
    #     "q_conc": 1e-9,
    #     # Final temperature of the annealing process.
    #     "temp": 21
    # }

    def __init__(self, **kwargs):

        # for arg, val in list(self.defaults.items()):
        #     setattr(self, arg, val)

        for arg, val in list(kwargs.items()):
            setattr(self, arg, val)

        # if isinstance(sess_or_client, cupyck.Client):
        #     self.client = sess_or_client
        #     self.session = None

        # elif isinstance(sess_or_client, Session):
        #     self.client = None
        #     self.session = sess_or_client

        # else:
        #     raise ValueError("must provide valid session or client")


    def simulate(self, feature_seq_pairs):
        """
        Takes a batch of pairs of feature sequences as a pandas dataframe, and
        estimates the CAS9 specificity for each

        """

#         with multiprocessing.Pool(8) as pool:
#             results = pool.map(f, feature_seq_pairs.iterrows())

#         return np.array(results)
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
            return np.max(site_scores, axis=1)
        else:
            return np.array(
                [crispr_specificity(p.target_features, p.query_features) for _, p in feature_seq_pairs.iterrows()]
            )