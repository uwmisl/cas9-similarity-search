import tensorflow as tf
from tensorflow.keras import layers

import numpy as np
import pandas as pd

class EncoderTrainer:
    """
    Glues together the models to ensure that similar images yield DNA sequences more likely to hybridize.
    Used to train the encoder such that it encodes these images accordingly.
    """

    def __init__(self, encoder, predictor):

        self.encoder = encoder
        self.predictor = predictor

        # Can't use the sequential Keras model anymore because we're combining two data streams (i.e. no longer strictly sequential).
        # Instead, we use the functional model.

        # Remember: Batch dimension is implied.
        # "2" - for pair of images.
        X_pairs = layers.Input([2, encoder.input_dim])
        # Split the images (since a Keras model can only take one input)
        # Slices: (batch dimension, first item in the pair, remaining feature vector dimensions)
        #  result -> (batch dimension, input dimensions)

        # Essentially, we started with a batch of feature-vector pairs...
        # ...And turned them into a pair of feature-vector batches.
        X1, X2 = layers.Lambda(lambda X: (X[:,0,:], X[:,1,:]))(X_pairs)

        distances = layers.Lambda(lambda Xs: tf.sqrt(tf.reduce_sum(tf.square(Xs[0]-Xs[1]), axis=1)))([X1,X2])

        # Independently transforms the batches of feature vectors into soft-max encoded DNA sequences.
        S1 = encoder(X1)
        S2 = encoder(X2)

        # Glue them back together! Back into a batch of feature vector pairs.
        S_pairs = layers.Lambda(
            lambda Ss: tf.stack(Ss, axis=-1)
        )([S1,S2])

        # Dimensions: (batch_size x 80 x 4 x 2 ) (i.e. batch size x DNA length x # of nucleotides x 2)

        # Swaps dimensions for the predictor, which wants (batch-size x 2 x DNA length x 4)
        S_pairs_T = layers.Lambda(lambda S: tf.transpose(S, [0, 3, 1, 2]))(S_pairs)

        # y_h: Estimated output
        y_h = layers.Lambda(tf.function(predictor))(S_pairs_T)
        y_h_T = layers.Reshape([1])(y_h)

        # Convenience model for internal inspection
        self.calcseq = tf.keras.Model(inputs=X_pairs, outputs=S_pairs_T)
        # Calcdists exists as a convenience property, if one needs to perform a distance calculation on the GPU at the same time (no training happening).
        self.calcdists = tf.keras.Model(inputs=X_pairs, outputs=distances)
        # The actual trainable model.
        self.model = tf.keras.Model(inputs=X_pairs, outputs=y_h_T)