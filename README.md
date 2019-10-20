# Coresets for dictionary learning

The following code implements the coreset for dictionary learning algorithm as it appears in

1. Feldman D., Feigin M. and Sochen N., [“Learning Big (Image) Data via Coresets for Dictionaries”](http://www.mit.edu/~michaf/publication/feldman-2013/JMIV2013.pdf), in Journal of Mathematical Imaging and Vision, March 2013
2. Feigin M., Feldman D. and Sochen N., [“From High Definition Image to Low Space Optimization”](http://www.mit.edu/~michaf/publication/feigin-2012-a/SSVM2012.pdf), in Scale Space and Variational Methods (SSVM), Ein-Gedi, Israel, 2011

The code for building a Coreset for KSVD dictionary learning is included in the following two files:

| | |
| --- | --- |
| KSVD Coreset code | `KSVDCoresetAlg.m` |
| SVD Coreset code | `SVDCoresetAlg.m` 
| | |

Note that the second file is only needed if a non constant first approximation is requested.

Two auxiliary functions that may be of use, decompose an image into a patches and build the vector matrix, and re-composition of the image. Note that these version do not create overlapping patches, which are required for some algorithms (mainly denoising) that use KSVD.

| | |
| --- | --- |
| Decompose image into patch vectors | `BuildY.m` |
| Rebuild image from vectors matrix | `YToImage.m` |
| | |

The code itself is easy to use.

To compute a dictionary using KSVD the KSVD-Box and OMP-Box code by Ron Rubinstein is needed. Original versions can be downloaded from [here](http://www.cs.technion.ac.il/~ronrubin/software.html). For reference see:

- Rubinstein, R., Zibulevsky, M., and Elad, M., ["Efficient implementation of the K-SVD algorithm using batch orthogonal matching pursuit"](http://cs.technion.ac.il/users/wwwb/cgi-bin/tr-get.cgi/2008/CS/CS-2008-08.revised.pdf) in CS Technion, 2008

Next, the easiest approach is to just compute the actual dictionary by:

1. Create a new KSVDCoresetAlg object

        K = KSVDCoresetAlg;

2. Change the default settings for K.sampleSize and K.svdVecs if needed
3. Decompose images into patches to get the input to the dictionary learning code (default is 8 X 8 patches)

        P = BuildY(img);

4. Compute the dictionary (note that the input vectors are given in the rows of P).

        D = K.computeDictionary(P, sparsity, dictsize, KSVD_iter_num);

5. Compute the reconstruction coefficients for the dictionary

        X = omp(D, P, D'*D, sparsity);

6. Reconstruct patches from dictionary and coefficients

        RP = D*X;

7. Reconstruct the image from the reconstructed patches (assumes default block size)

        P = YToImage(RP, size(img));

The Coreset may be obtained directly by using the K.computeCoreset function, but note the implementation in K.computeDictionary as to the correct usage of the approximated dictionary (variable tD in K.computeDictionary)
