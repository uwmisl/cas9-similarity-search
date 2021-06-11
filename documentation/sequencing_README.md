sequencing_index,encoder,query_id,qt_ratio,dsDNA_conc,n_randomers,rep

# Sequencing_index

The sequencing index (e.g. C2, C5, C6) refers to a unique sample.

# Encoder

TODO: callie q: What's the encoder in this context? Does this uniquely identify the encoder trained at a given date?  https://github.com/uwmisl/cas9-similarity-search/issues/2

# Query_ID

TODO: callie q: check for correctness?  https://github.com/uwmisl/cas9-similarity-search/issues/2
The label for the query used in the particular experiment. This is the file name for the query image.

# QT Ratio

The query:target ratio. This ratio is determined by the number of query molecules: target molecules.
A subtlety: the number of target molecules is quantified at the dsDNA step (right before the 3 PCR cycles of linear amplification).


# dsDNA_conc

TODO: callie q: check for correctness?  https://github.com/uwmisl/cas9-similarity-search/issues/2
The double stranded DNA concentration in nM. The dsDNA is the partially double stranded target immediately before the query is bound to it.

# n_randomers

TODO: callie q: check for correctness?  https://github.com/uwmisl/cas9-similarity-search/issues/2
The number of the randomers added to the experiment. Each randomer has the normal conserved regions and the 80nt feature region is randomized. The 30nt ID region is also randomized.

# rep

The replicate number. I.e., the unique identifier for the number of times this exact experiment has been done.