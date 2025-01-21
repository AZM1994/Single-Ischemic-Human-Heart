## General Description 

Mutations that occur during cancer development can be classified into different mutational processes. This is done using the mutations’ sequence contexts and statistical tools that decompose the observed mutation spectra into statistically independent mutational components, or "signatures". Many so-called refitting algorithms have been developed to find a linear combination of these independent signatures that can describe the data and with this we can assign mutations from a given tumor to the underlying mutational processes. However, these algorithms were developed years ago, when a set of 30 signatures was originally identified, and they have not been tested for the newest catalogs that contain more than 70 different signatures. Furthermore, estimated per-tumor signature weights are noisy, yet their errors are rarely quantified in existing algorithms. 

We have developed a method that uses Deep Neural Networks to do signature refitting. This algorithm, called SigNet, is easily adaptable to any new catalog of mutational signatures and provides prediction intervals for the contribution of each signature in a single tumor. This showcases how, using Artificial Neural Networks, we can take advantage of the correlations between the different mutational processes that occur during carcinogenesis to obtain accurate signature decompositions, even when the number of mutations in the sampled tumor is very low. For more details on how the method works please see [link to the paper].

### What is signature refitting?

Given a set of known mutational signatures, signature refitting methods aim to find a linear combination of these signatures that can reconstruct the mutational profile of a given sample. The mutational signatures have been identified using the conventional 96 mutation type classification, considering not only the mutated base (six substitution subtypes: C>A, C>G, C>T, T>A, T>C, and T>G), but also the bases immediately 5’ and 3’. The currently most used set of mutational signatures is the one provided by COSMIC, [[1]](#1). SigNet uses COSMIC v3.1, which contains 72 signatures. 

## References
<a id="1">[1]</a> 
COSMIC. Mutational Signatures (v3.3 - June 2022).
Accessed on 14th September 2023.
https://cancer.sanger.ac.uk/signatures/sbs/
