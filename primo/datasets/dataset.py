import abc
import sklearn.metrics
import numpy as np

class Dataset(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def random_pairs(self, batch_size):
        pass

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
            distance_matrix = sklearn.metrics.pairwise_distances(anchor_vals[search_rows], chunk_vals)
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