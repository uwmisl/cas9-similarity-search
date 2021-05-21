import pandas as pd

import tensorflow as tf
from tensorflow.keras import layers

from . import _default_sequences
from ..tools import sequences as seqtools

def entropy_regularizer(strength):
    # TODO: Verify with callie this is correct: https://github.com/uwmisl/cas9-similarity-search/issues/2
    """
    This function adds some regularization to the model, meaning it penalizes the model for
    getting 'too complex' so that it learns to chill out. How harshly you punish the model for
    complexity is mediated by the 'strength' parameter.
    Parameters
    ----------
    strength : float
        How drastically the model should correct itself when regularizing.
    """

    def encoder_entropy(seq_probs):
        seq_probs.shape.assert_is_compatible_with([None, None, 4])

        ent_by_position = -tf.reduce_sum(
            seq_probs * tf.log(seq_probs + 1e-10),
            axis = 2
        )
        mean_ent_by_sequence = tf.reduce_mean(
            ent_by_position,
            axis = 1
        )
        mean_ent_by_batch = tf.reduce_mean(
            mean_ent_by_sequence,
            axis = 0
        )

        return strength * mean_ent_by_batch

    return encoder_entropy

class Encoder:

    defaults = {
        # TODO: Verify with callie this is correct: https://github.com/uwmisl/cas9-similarity-search/issues/2
        # The reason the default encoder input is a 4096-dimensional vector is
        # because we're representing our images through an embedding that was learned
        # by a computer vision model known as VGG [1]. We're borrowing the output of the
        # 2nd fully-connected layer (i.e. the FC2), which spits out a 4096-by-1 vector.
        #
        # If you're very curious about VGG's innerworkings, you can see an example tensorflow
        # implementation here [2, 3].
        #
        # Note for future users: If you ever decide to use a different model VGG16,
        # you'd probably want to change the input dimension here.
        #
        # [1] - https://neurohive.io/en/popular-networks/vgg16/
        # [2] - https://www.cs.toronto.edu/~frossard/post/vgg16/
        # [3] - https://github.com/kentsommer/VGG16-Image-Retrieval/blob/master/vgg16_example.py#L237
        #
        "input_dim": 4096,

        # TODO: Verify with callie this is correct: https://github.com/uwmisl/cas9-similarity-search/issues/2
        # The feature region of our engineered DNA sequence is 80 nucleotides long.
        "output_len": 80,

        # Regularization penalty post softmax and helps prevent overfitting.
        # This value, 1e-2, was experimentally determined.
        # Since the output was softmax, a valid range of regularization strength is between 0 and 1.
        "entropy_reg_strength": 1e-2
    }

    def __init__(self, model_path = None, **kwargs):
        """DNA Sequence Encoder. This class handles

        TODO: Verify with callie that our documentation assumptions
        are correct. https://github.com/uwmisl/cas9-similarity-search/issues/2

        Parameters
        ----------
        model_path : [type], optional
            [description], by default None
        """
        for arg, val in _default_sequences.items():
            setattr(self, arg, val)

        for arg, val in self.defaults.items():
            setattr(self, arg, val)

        for arg, val in kwargs.items():
            setattr(self, arg, val)

        # TODO: Verify with callie that this '4' is the number of DNA bases. https://github.com/uwmisl/cas9-similarity-search/issues/2
        number_of_bases = len(seqtools.bases)
        if model_path is None:
            self.model = tf.keras.Sequential([
                layers.Dense(self.input_dim/2, activation = 'relu', input_shape=[self.input_dim]),
                layers.Dense(self.output_len * number_of_bases, activation='relu'),
                layers.Reshape([self.output_len, number_of_bases]),
                layers.Activation('softmax'),
                layers.Lambda(

                    # TODO: Verify with callie this is correct: https://github.com/uwmisl/cas9-similarity-search/issues/2
                    # Just using the identity because we don't want to transform the softmaxxed output,
                    # we just want to make sure we learn an output encoding that's regularized (i.e. not crazy complex/over-fitting)
                    lambda x: x,

                    # TODO: Verify with callie this is correct: https://github.com/uwmisl/cas9-similarity-search/issues/2
                    # Using an "entropy" regulator because we passed the output through a softmax.
                    activity_regularizer=entropy_regularizer(
                        self.entropy_reg_strength
                    )
                )
            ], name='encoder')

        else:
            self.model = tf.keras.models.load_model(model_path)

    def __call__(self, X):
        return self.model(X)

    def trainable(self, is_trainable):
        # TODO: Python 3: Perhaps we could make this a property? https://github.com/uwmisl/cas9-similarity-search/issues/4
        # TODO: Verify with callie this is correct: https://github.com/uwmisl/cas9-similarity-search/issues/2
        """Whether to use this model in training mode or inference mode.

        Parameters
        ----------
        is_trainable : boolean
            Whether this model should be used for training or inference.
        """
        for layer in self.model.layers:
            layer.trainable = is_trainable

    def encode_feature_seqs(self, X):
        onehots = self.model.predict(X)
        return seqtools.onehots_to_seqs(onehots)

    # todo
    def build_targets(self, feature_seqs, barcoder):
        pass

    # todo
    def build_queries(self, feature_seqs):
        pass


    # evaluate ???

    def save(self, model_path):
        self.model.save(model_path)

