sequencing_index,encoder,query_id,qt_ratio,dsDNA_conc,n_randomers,rep

# Sequencing_index

The sequencing index (e.g. C2, C5, C6) refers to a unique sample.

# Encoder

Uniquely identifies the encoder trained at a given date.

# Query_ID

The label for the query used in the particular experiment. This is the file name for the query image.

# QT Ratio

The query:target ratio. This ratio is determined by the number of query molecules: target molecules.
A subtlety: the number of target molecules is quantified at the dsDNA step (right before the 3 PCR cycles of linear amplification).


# dsDNA_conc

TODO: yuan q: check for correctness?  https://github.com/uwmisl/cas9-similarity-search/issues/9
The double stranded DNA concentration in nM. The dsDNA is the partially double stranded target immediately before the query is bound to it.

# n_randomers

TODO: yuan q: check for correctness? https://github.com/uwmisl/cas9-similarity-search/issues/9

The number of the randomers added to the experiment. Each randomer has the normal conserved regions and the 80nt feature region is randomized. The 30nt ID region is also randomized.

# rep

TODO: yuan q: check for correctness? https://github.com/uwmisl/cas9-similarity-search/issues/9

The replicate number. I.e., the unique identifier for the number of times this exact experiment has been done.