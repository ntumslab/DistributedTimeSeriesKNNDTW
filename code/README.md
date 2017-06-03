# Source Code

## Compilation (GCC for example)

```
gcc -o code -O3 -std=c99 -lm code.c
```

## Running

```
./code arguments -f input_file
```

## Arguments

* -f: dataset file name (Required)
* -k: the number of candidates to be retrieved (Default 10)
  * 1 <= k <= s
* -s: the number of sensors (candidates) (Default 500)
  * 1 <= s <= Total number of time series in a dataset - 1
  * Able to test multiple values (e.g. 500,1000,1500) split by comma without space
* -m: the number of sites (Default 50)
  * 1 <= m <= s
  * Able to test multiple values (e.g. 50,100,150) split by comma without space
* -t: time series length (Default maximum length)
  * 1 <= t <= Maximum length provided by a dataset
* -q: sensor (time series) index to be assigned as the query (Default -1)
  * -2 <= q <= s
  * q == -1: Let every time series be the query in turns (for small data)
  * q == -2: Randomly choose a time series as the query for 100 times (for big data)
* -r: The parameter R. The number of segments cut from a segment at the previous level (Default 2)
  * 2 <= r
* -n: The parameter N. The number of segments at the first level (Default round(sqrt(t / 2)))
  * 1 <= n <= t

## Input file

The format of **input_file** follows the setting of UCR datasets:

```
Keogh, E., Zhu, Q., Hu, B., Hao. Y.,  Xi, X., Wei, L. & Ratanamahatana, C. A.
(2011). The UCR Time Series Classification/Clustering Homepage:
http://www.cs.ucr.edu/~eamonn/time_series_data/
```

## Output file

The file extension name of **input_file** is replaced with "out".
The file contains the following group of tables.

Table format:
* **Row**: number of sensors (candidates)
* **Column**: number of sites
* **Entry**: bandwidth ratio = framework bandwidth (bits) / naive approach bnadwidth (bits)

Table contents:
* **Table 1**: mean of bandwidth ratio
* **Table 2**: standard deviation of bandwidth ratio
* **Table 3**: mean of candidate time series used for initialization
* **Table 4**: mean of candidate time series pruned by the level-1 query in the server
* **Table 5**: mean of candidate time series pruned by LB_KimFL
* **Table 6**: mean of candidate time series pruned by LB_Keogh
* **Table 7**: mean of candidate time series pruned by LB_MS
* **Table 8**: mean of candidate time series pruned unable to be pruned
* **Table 9**: mean of compression ratio of the baseline model (if COMPRESSION macro name is defined, or 0 if not).
* **Table 10**: mean of compression ratio of the framework model (if COMPRESSION macro name is defined, or 0 if not).
* **Table 11**: mean of candidate sites used for initialization
* **Table T = 12 to 20**: mean of candidate sites pruned at the resolution level (T - 11)
* **Table 21**: mean of candidate sites pruned at the resolution level 10 or more

## Experimented Parameters for UCR Datasets

```
./code -s 500,1000,1500 -m 50,100,150 -f [UCR dataset text file]
```


## Experimented Parameters for the Synthetic Dataset

```
./code -s 12499 -m 500,1000,1500 -k 30 -q -2 -f synthetic_data.txt
```

## Other Experiment Scenarios

* **Equal-size segmentation**: cancel the comment to the macro "EQUAL_SIZE" at the head of code.c
* **Oracle initialization**: cancel the comment to the macro "ORACLE_INITIALIZATION" at the head of code.c
* **Random initialization**: cancel the comment to the macro "RANDOM_INITIALIZATION" at the head of code.c
* **Compression**: cancel the comment to the macro "COMPRESSION" at the head of code.c, but the OS where code.c runs should support the "zip" command.
