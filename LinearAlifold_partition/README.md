
# LinearAlifold: Linear-Time Consensus Structure Prediction for RNA Alignments

Liang Zhang, Sizhen Li, He Zhang, David H. Mathews, Liang Huang*

\* corresponding author


## Dependencies
gcc 4.8.5 or above; 
python2.7

## To Compile
```
make
```

## To Run
(input: a Multiple Sequence Alignment (MSA)):
```
cat MSA_file | ./linearalifold [OPTIONS]
```

OPTIONS:
```
-b BEAM_SIZE
```
The beam size (default 100). Use 0 for infinite beam.
```
-o FILE_NAME
```
Outputs base pairing probability matrix to a file with user specified name. (default False)
```
-r FILE_NAME
```
Output base pairing probability matrix to a file with user specified name (overwrite if the file exists). (default False)
```
--prefix PREFIX_NAME
```
Outputs base pairing probability matrices to files with user specified prefix. (default False)
```
-p
```
Partition function calculation only. (default False)
```
-c CUTOFF
```
Only output base pair probability larger than user specified threshold between 0 and 1. (DEFAULT=0.0)
```
--dumpforest
```
dump forest (all nodes with inside [and outside] log partition functions but no hyperedges) for downstream tasks such as sampling and accessibility (DEFAULT=None)

```
--mea or -M
```
get MEA structure, (DEFAULT=FALSE)

```
--gamma
```
set MEA gamma, (DEFAULT=3.0)

```
--bpseq
```
output MEA structure(s) in bpseq format instead of dot-bracket format

```
--mea_prefix
```
output MEA structure(s) to file(s) with user specified prefix name

```
--threshknot or -T
```
get ThreshKnot structure, (DEFAULT=FALSE)

```
--threshold
```
set ThreshKnot threshknot, (DEFAULT=0.3)

```
--threshknot_prefix
```
output ThreshKnot structure(s) to file(s) with user specified prefix name (default False)


## Example: Run Predict
```
cat alignment_fasta.fa | ./linearalifold
Free Energy of Ensemble: -26.43 kcal/mol
```

## Example: Run Partition Function Calculation Only
```
cat alignment_fasta.fa | ./linearalifold -p --verbose
beam size: 100
Free Energy of Ensemble: -26.43 kcal/mol
Partition Function Calculation Time: 0.00 seconds.
```

## Example: Run Prediction and Output MEA structure
```
cat alignment_fasta.fa | ./linearalifold -M
Free Energy of Ensemble: -26.43 kcal/mol
...((.(((((((((....((...))))))))))).))
```

## Example: Run Prediction and Output ThreshKnot structure in bpseq format
The neotides in the output are from the first aligned sequence.
```
cat alignment_fasta.fa | ./linearalifold -T --threshold 0
Free Energy of Ensemble: -26.43 kcal/mol
1 C 0
2 U 0
3 C 0
4 A 38
5 C 37
6 A 0
7 A 35
8 C 34
9 G 33
10 U 32
11 U 31
12 U 30
13 G 29
14 U 28
15 G 27
16 C 0
17 C 0
18 U 0
19 C 0
20 A 0
21 G 26
22 U 0
23 U 36
24 A 0
25 C 0
26 C 21
27 C 15
28 G 14
29 U 13
30 A 12
31 G 11
32 A 10
33 U 9
34 G 8
35 U 7
36 A 23
37 G 5
38 U 4
```
