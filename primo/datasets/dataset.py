import abc
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



# TODO: I moved the `keras_batch_generator` function to `primo/datasets/dataset.py` based on the principle of Cohesion-- because it feels like this code is closely tied to the abstract dataset generators. callie, if there was a specific reason for keeping this function in this notebook in this cell, please correct me. https://github.com/uwmisl/cas9-similarity-search/issues/2
def keras_batch_generator(dataset_batch_generator, similarity_threshold):
    # Yield datasets
    # TODO: Verify with callie this understanding is correct https://github.com/uwmisl/cas9-similarity-search/issues/2
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