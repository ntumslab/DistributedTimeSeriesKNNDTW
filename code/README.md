k nearest neighbor search on distributed time series

====

Compilation (GCC for example)

gcc -o code -O3 -std=c99 -lm code.c

====

Running

./code [parameters] -f [file_name]

====

Input file: [file_name]
The file format follows the setting of UCR datasets:

Keogh, E., Zhu, Q., Hu, B., Hao. Y.,  Xi, X., Wei, L. & Ratanamahatana, C. A.
(2011). The UCR Time Series Classification/Clustering Homepage:
www.cs.ucr.edu/~eamonn/time_series_data/

====

Output file: [file_name].out, the original file extension will be removed if [file_name] has it.
The file is a set of tables.

Row: Number of sensors (candidates)
Column: Number of sites
Cell: Bandwidth ratio = framework bandwidth (bits) / naive approach bnadwidth (bits)

Table 1: Mean of bandwidth ratio
Table 2: Standard deviation of bandwidth ratio
Table 3: Mean of candidate time series used for initialization
Table 4: Mean of candidate time series pruned by the level-1 query in the server
Table 5: Mean of candidate time series pruned by LB_KimFL
Table 6: Mean of candidate time series pruned by LB_Keogh
Table 7: Mean of candidate time series pruned by LB_MS
Table 8: Mean of candidate time series pruned unable to be pruned
Table 9: Mean of compression ratio of the baseline model (if COMPRESSION macro name is defined, or 0 if not).
Table 10: Mean of compression ratio of the framework model (if COMPRESSION macro name is defined, or 0 if not).
Table 11: Mean of candidate sites used for initialization
Table T = 12 to 20: Mean of candidate sites pruned at the resolution level (T - 11)
Table 21: Mean of candidate sites pruned at the resolution level 10 or more

====

Parameters

-f: Dataset file name (Required)

-k: The number of candidates to be retrieved (Default 10)
	1 <= k <= s

-s: The number of sensors (candidates) (Default 500)
	1 <= s <= Total number of time series in a dataset - 1
	Able to test multiple values (e.g. 500,1000,1500) split by comma without space

-m: The number of sites (Default 50)
	1 <= m <= s
	Able to test multiple values (e.g. 50,100,150) split by comma without space

-t: Time series length (Default maximum length)
	1 <= t <= Maximum length provided by a dataset

-q: Sensor (time series) index to be assigned as the query (Default -1)
	-2 <= q <= s
	q == -1: Let every time series be the query in turns (for small data)
	q == -2: Randomly choose a time series as the query for 100 times (for big data)

-r: The parameter R. The number of segments cut from a segment at the previous level (Default 2)
	2 <= r

-n: The parameter N. The number of segments at the first level (Default round(sqrt(t / 2)))
	1 <= n <= t

====

Experiments for the UCR datasets (Most parameters use the default settings)

./code -s 500,1000,1500 -m 50,100,150 -f [UCR dataset text file]

====

Experiments for the synthetic dataset

./code -s 12499 -m 500,1000,1500 -k 30 -q -2 -f synthetic_data.txt

====

Other scenarios

Equal-size segmentation: Cancel the comment to the macro "EQUAL_SIZE" at the head of code.c
Oracle initialization: Cancel the comment to the macro "ORACLE_INITIALIZATION" at the head of code.c
Random initialization: Cancel the comment to the macro "RANDOM_INITIALIZATION" at the head of code.c
Compression: Cancel the comment to the macro "COMPRESSION" at the head of code.c, but the OS where code.c runs should support the "zip" command.
