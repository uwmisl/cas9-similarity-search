# Utilities for location the file paths of import files.

from pathlib import Path

def get_encoder_model_path(isDocker: bool = False) -> str:
    """Gets the filepath of a stored DNA encoder.

    Parameters
    ----------
    isDocker : bool
        Whether this program is being run from a Docker container (and thus has different file paths).

    Returns
    -------
    str
        Absolute path to the DNA encoder.

    Raises
    ------
    NotImplementedError
        Non-Docker file paths haven't been implemented yet.
    """
    path = "/tf/primo/data/models/encoder-model.h5"
    if not isDocker:
        # TODO: Replace with a method for returning an on-disk filepath: https://github.com/uwmisl/cas9-similarity-search/issues/6
        raise NotImplementedError("TODO: Replace with a method for returning filepath https://github.com/uwmisl/cas9-similarity-search/issues/6")
    return path

def get_query_features_path(isDocker: bool = False) -> str:
    """Gets the filepath for the set of query features.

    Parameters
    ----------
    isDocker : bool
        Whether this program is being run from a Docker container (and thus has different file paths).

    Returns
    -------
    str
        Absolute path to the query features.

    Raises
    ------
    NotImplementedError
        Non-Docker file paths haven't been implemented yet.
    """
    path = "/tf/primo/data/queries/features.h5"
    if not isDocker:
        # TODO: Replace with a method for returning an on-disk filepath: https://github.com/uwmisl/cas9-similarity-search/issues/6
        raise NotImplementedError("TODO: Replace with a method for returning filepath https://github.com/uwmisl/cas9-similarity-search/issues/6")
    return path


def get_encoded_query_sequences_path(isDocker: bool = False) -> str:
    """Where to save a set of encoded query sequences.

    Parameters
    ----------
    isDocker : bool
        Whether this program is being run from a Docker container (and thus has different file paths).

    Returns
    -------
    str
        Absolute path to the encoded query sequences features.

    Raises
    ------
    NotImplementedError
        Non-Docker file paths haven't been implemented yet.
    """
    path = "/tf/primo/data/queries/feature_seqs.h5"
    if not isDocker:
        # TODO: Replace with a method for returning an on-disk filepath: https://github.com/uwmisl/cas9-similarity-search/issues/6
        raise NotImplementedError("TODO: Replace with a method for returning filepath https://github.com/uwmisl/cas9-similarity-search/issues/6")
    return path


def get_distance_store_path(isDocker: bool = False) -> str:
    """Where to save the cache of Euclidean distances between queries and targets.

    Parameters
    ----------
    isDocker : bool
        Whether this program is being run from a Docker container (and thus has different file paths).

    Returns
    -------
    str
        Absolute path to the target-query feature distances.

    Raises
    ------
    NotImplementedError
        Non-Docker file paths haven't been implemented yet.
    """
    path = "/tf/primo/data/targets/query_target_dists.h5"
    if not isDocker:
        # TODO: Replace with a method for returning an on-disk filepath: https://github.com/uwmisl/cas9-similarity-search/issues/6
        raise NotImplementedError("TODO: Replace with a method for returning filepath https://github.com/uwmisl/cas9-similarity-search/issues/6")
    return path

def get_sequence_store_path(isDocker: bool = False) -> str:
    """Where to save the DNA sequence encoding of target features.

    Parameters
    ----------
    isDocker : bool
        Whether this program is being run from a Docker container (and thus has different file paths).

    Returns
    -------
    str
        Absolute path to the  DNA sequence encoding of target features.

    Raises
    ------
    NotImplementedError
        Non-Docker file paths haven't been implemented yet.
    """
    path = "/tf/primo/data/targets/feature_seqs.h5"
    if not isDocker:
        # TODO: Replace with a method for returning an on-disk filepath: https://github.com/uwmisl/cas9-similarity-search/issues/6
        raise NotImplementedError("TODO: Replace with a method for returning filepath https://github.com/uwmisl/cas9-similarity-search/issues/6")
    return path

def get_target_feature_path(prefix: str, isDocker: bool = False) -> str:
    """Path to a set of target features, with a prefix.

    Parameters
    ----------
    prefix : str
        Identifier for target feature set. Often times the target images are too large to store in a single file, and thus are split across
        multiple ones, each file being given some unique identifier (e.g. hexadecimal index)

    isDocker : bool
        Whether this program is being run from a Docker container (and thus has different file paths).

    Returns
    -------
    str
        Absolute path to a particular target feature set.

    Raises
    ------
    NotImplementedError
        Non-Docker file paths haven't been implemented yet.
    """
    path = f"/tf/open_images/targets/features/targets_{prefix}.h5"
    if not isDocker:
        # TODO: Replace with a method for returning an on-disk filepath: https://github.com/uwmisl/cas9-similarity-search/issues/6
        raise NotImplementedError("TODO: Replace with a method for returning filepath https://github.com/uwmisl/cas9-similarity-search/issues/6")
    return path

def get_extended_target_distance_store_path(isDocker: bool = False) -> str:
    """Where to save the query-target distances in the extended data set.
    'Extended' in this context refers to the ~9 million images not used in either the 'target' or the 'train' set.

    Parameters
    ----------
    isDocker : bool
        Whether this program is being run from a Docker container (and thus has different file paths).

    Returns
    -------
    str
        Absolute path to the query-target distances for the extended targets.

    Raises
    ------
    NotImplementedError
        Non-Docker file paths haven't been implemented yet.
    """
    path = "/tf/primo/data/extended_targets/query_target_dists.h5"
    if not isDocker:
        # TODO: Replace with a method for returning an on-disk filepath: https://github.com/uwmisl/cas9-similarity-search/issues/6
        raise NotImplementedError("TODO: Replace with a method for returning filepath https://github.com/uwmisl/cas9-similarity-search/issues/6")
    return path

def get_extended_sequence_store_path(isDocker: bool = False) -> str:
    """Where to cache the encoded DNA sequences for the extended target features.
    'Extended' in this context refers to the ~9 million images not used in either the 'target' or the 'train' set.

    Parameters
    ----------
    isDocker : bool
        Whether this program is being run from a Docker container (and thus has different file paths).

    Returns
    -------
    str
        Absolute path to the DNA sequence encoding of the extended target features.

    Raises
    ------
    NotImplementedError
        Non-Docker file paths haven't been implemented yet.
    """
    path = "/tf/primo/data/extended_targets/feature_seqs.h5"
    if not isDocker:
        # TODO: Replace with a method for returning an on-disk filepath: https://github.com/uwmisl/cas9-similarity-search/issues/6
        raise NotImplementedError("TODO: Replace with a method for returning filepath https://github.com/uwmisl/cas9-similarity-search/issues/6")
    return path

def get_extended_target_feature_path(prefix: str, isDocker: bool = False) -> str:
    """Path to a set of extended target features, with a prefix.
    'Extended' in this context refers to the ~9 million images not used in either the 'target' or the 'train' set.

    Parameters
    ----------
    prefix : str
        Identifier for extended target feature set. Often times the target images are too large to store in a single file, and thus are split across
        multiple ones, each file being given some unique identifier (e.g. hexadecimal index)

    isDocker : bool
        Whether this program is being run from a Docker container (and thus has different file paths).

    Returns
    -------
    str
        Absolute path to a particular extended target feature set.

    Raises
    ------
    NotImplementedError
        Non-Docker file paths haven't been implemented yet.
    """
    path = f"/tf/open_images/extended_targets/features/extended_targets_{prefix}.h5"
    if not isDocker:
        # TODO: Replace with a method for returning an on-disk filepath: https://github.com/uwmisl/cas9-similarity-search/issues/6
        raise NotImplementedError("TODO: Replace with a method for returning filepath https://github.com/uwmisl/cas9-similarity-search/issues/6")
    return path