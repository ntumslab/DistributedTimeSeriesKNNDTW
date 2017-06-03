# Distributed Time Series KNN DTW

Supplementary material and source code of paper "Bandwidth-Efficient Distributed k-Nearest-Neighbor Search with Dynamic Time Warping"

## Abstract

```
We study the fundamental k-nearest neighbor (kNN) search problem on distributed time series. A server has constantly received various reference time series Q of length X and seeks the exact kNN over a collection of time series distributed across a set of M local sites. When X and M are large, and when the amount of query increases, simply sending each Q to all M sites incurs high communication bandwidth costs, which we would like to avoid. Prior work has presented a communication-efficient kNN algorithm for the Euclidean distance similarity measure. In this paper, we present the first communication-efficient kNN algorithm for the dynamic time warping (DTW) similarity measure, which is generally believed a better measure for time series. To handle the complexities of DTW, we design a new multi-resolution structure for the reference time series, and multi-resolution lower bounds that can effectively prune the search space. We present a new protocol between the server and the local sites that leverages multi-resolution pruning for communication efficiency and cascading lower bounds for computational efficiency. Empirical studies on both real-world and synthetic data sets show that our method reduces communication bandwidth by up to 92%.
```

## References

```
Hsu, Chin-Chi, et al. "Bandwidth-efficient distributed k-nearest-neighbor search with dynamic time warping." Big Data (Big Data), 2015 IEEE International Conference on. IEEE, 2015.
```

IEEE link (to download BibTex reference format):
http://ieeexplore.ieee.org/document/7363799/?reload=true&arnumber=7363799

## Directories

* **supplementary_materials**: additional proofs or explanations about our proposed framework.
* **code**: source code for our experiments.
* **data**: datasets for our experiments.
