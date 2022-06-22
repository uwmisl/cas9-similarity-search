# Files in this directory are meant to be run in (or are the product from) the [`primo-similarity-search`](https://github.com/uwmisl/primo-similarity-search) github repository.

## What's the purpose of this directory?
This directory exists to compare the ability of DNA hybridization to simulate the performance of similarity search when there is an 80nt space (as presented in Bee et al. 2021) or a 20nt space (generated using all the same tools as Bee et al.).

## What's in this directory?
`simulate_pairs_hyb.ipynb` is the notebook run in the `primo-similarity-search` repository that generated the comparison analysis. It requires the use of the `.h5` files included here (generated in the `primo-similarity-search` repo by only changing the length from 80 to 20).


