import abc
# import sklearn.metrics
import numpy as np
import tensorflow as tf
import os
import sys
import pandas as pd

def pairwise_dist (A, B):
  """
  Computes pairwise distances between each elements of A and each elements of B.
  Args:
    A,    [m,d] matrix
    B,    [n,d] matrix
  Returns:
    D,    [m,n] matrix of pairwise distances
  """
  with tf.compat.v1.variable_scope('pairwise_dist'):
    # squared norms of each row in A and B
    na = tf.reduce_sum(tf.square(A), 1)
    nb = tf.reduce_sum(tf.square(B), 1)

    # na as a row and nb as a co"lumn vectors
    na = tf.reshape(na, [-1, 1])
    nb = tf.reshape(nb, [1, -1])

    # return pairwise euclidead difference matrix
    D = tf.sqrt(tf.maximum(na - 2*tf.matmul(A, B, False, True) + nb, 0.0))
  return D

class Dataset(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def random_pairs(self, batch_size):
        pass

    def random_features(self, batch_size):
        feature_dir = os.path.join(self.path, 'features')
        files = os.listdir(feature_dir)

        while True:

            f_a, f_b = np.random.choice(files, 2, replace=False)
            sys.stdout.write("switching to %s and %s\n" % (f_a, f_b))

            df1 = pd.read_hdf(os.path.join(feature_dir, f_a))
            df2 = pd.read_hdf(os.path.join(feature_dir, f_b))

            df = pd.concat([df1, df2])
            n = len(df)

            for _ in range(self.switch_every):
                pairs = np.random.permutation(n)[:batch_size]
                yield df.index.values[pairs], df.values[pairs]

    def balanced_pairs(self, batch_size, sim_thresh):

        pair_generator = self.random_pairs(batch_size)

        while True:

            n_batch = 0
            batch_ids  = []
            batch_vals = []

            while n_batch < batch_size:
                chunk_ids, chunk_vals = next(pair_generator)

                distances = np.sqrt(
                    np.square(chunk_vals[:,0] - chunk_vals[:,1]).sum(1)
                )

                similar = distances <= sim_thresh
                n_sim = similar.sum()

                batch_ids.extend([
                    chunk_ids[similar],
                    chunk_ids[~similar][:n_sim]
                ])

                batch_vals.extend([
                    chunk_vals[similar],
                    chunk_vals[~similar][:n_sim]
                ])

                n_batch += 2 * n_sim

            batch_ids = np.concatenate(batch_ids)
            batch_vals = np.concatenate(batch_vals)

            perm = np.random.permutation(len(batch_vals))[:batch_size]

            yield batch_ids[perm], batch_vals[perm]


class Static(Dataset):

    def __init__(self, X):
        self.X = X

    def random_pairs(self, batch_size):
        n,d = self.X.shape
        while True:
            pairs = np.random.permutation(n)[:batch_size*2].reshape(-1,2)
            yield pairs, self.X[pairs]


def triplet_batch_generator(dataset_batch_generator, similarity_threshold):
    while True:
        n_batch = 0

        # Get a batch worth of anchors
        anchor_ids, anchor_vals = next(dataset_batch_generator)


        pos_filled = np.zeros_like(anchor_ids, dtype=bool)
        pos_ids = np.zeros_like(anchor_ids)
        pos_vals = np.zeros_like(anchor_vals)
        neg_filled = np.zeros_like(anchor_ids, dtype=bool)
        neg_ids = np.zeros_like(anchor_ids)
        neg_vals = np.zeros_like(anchor_vals)

        batches_searched = 0
        # Keep pulling batches and taking random positive and negative samples until all anchors both
        while (not np.all(pos_filled) or not np.all(neg_filled)) and batches_searched < 100:
            chunk_ids, chunk_vals = next(dataset_batch_generator)

            batches_searched += 1

            # Get index of samples which do not yet have a pos and negative pairing
            search_rows = np.invert(pos_filled) | np.invert(neg_filled)

            # Compute distance for those rows
            #distance_matrix = sklearn.metrics.pairwise_distances(anchor_vals[search_rows], chunk_vals)
            distance_matrix = pairwise_dist(
                tf.convert_to_tensor(anchor_vals[search_rows], dtype=np.float32),
                tf.convert_to_tensor(chunk_vals, dtype=np.float32),
            ).numpy()
            similar_matrix = distance_matrix <= similarity_threshold

            # Create a matrix of random numbers, used to randomly select among multiple matches for each anchor
            random_matrix = np.random.randint(1, similar_matrix.shape[1], size=similar_matrix.shape)

            # Select positive samples
            maxidx = np.argmax(similar_matrix * random_matrix, axis=1)
            maxval = np.max(similar_matrix * random_matrix, axis=1)
            update_idx = (maxval > 0)
            pos_filled[np.flatnonzero(search_rows)[update_idx]] = True
            pos_ids[np.flatnonzero(search_rows)[update_idx]] = chunk_ids[maxidx[update_idx]]
            pos_vals[np.flatnonzero(search_rows)[update_idx]] = chunk_vals[maxidx[update_idx]]

            # Select negative samples
            disimilar_matrix = ~similar_matrix
            maxidx = np.argmax(disimilar_matrix * random_matrix, axis=1)
            maxval = np.max(disimilar_matrix * random_matrix, axis=1)
            update_idx = (maxval > 0)
            neg_filled[np.flatnonzero(search_rows)[update_idx]] = True
            neg_ids[np.flatnonzero(search_rows)[update_idx]] = chunk_ids[maxidx[update_idx]]
            neg_vals[np.flatnonzero(search_rows)[update_idx]] = chunk_vals[maxidx[update_idx]]


        # It's not usually possible to find positive pairings for all anchors.
        # In this cases, use the anchor as its own positive pair
        pos_empty = np.invert(pos_filled)
        pos_vals[pos_empty] = anchor_vals[pos_empty]
        pos_ids[pos_empty] = anchor_ids[pos_empty]

        ids = np.stack([anchor_ids, pos_ids, neg_ids], axis=1)
        triplets = np.stack([anchor_vals, pos_vals, neg_vals], axis=1)

        yield ids, triplets

def keras_batch_generator(dataset_batch_generator, similarity_threshold):
    # Yield datasets
    while True:
        # This tuple contains:
        # indices: a positive integer uniquely identifying an image. This index is obtained by enumerating all the images in the dataset (before splitting them into test/train/validate datasets)
        # pairs:
        indices, pairs = next(dataset_batch_generator)
        # The Euclidean distances between the two vectors in each pair
        distances = np.sqrt(np.square(pairs[:,0,:] - pairs[:,1,:]).sum(1))
        # Whether or not the images in this pair should be considered 'similar'. This is a boolean value, represented by an int (0 or 1), and is determined by whether the aforementioned Euclidean distances between image feature vectors are under some pre-deterined "similarity threshold".
        similar = (distances < similarity_threshold).astype(int)
        # Yield a pair of sequences, and 0-or-1 indicating whether they're similar.
        yield pairs, similar