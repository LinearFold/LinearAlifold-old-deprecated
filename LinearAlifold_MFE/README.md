
# LinearAlifold: Linear-Time Consensus Structure Prediction for RNA Alignments

Liang Zhang, Sizhen Li, He Zhang, David H. Mathews, Liang Huang*

\* corresponding author


## Dependencies
GCC 4.8.5 or above; 
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


--verbose
```
Prints out runtime in seconds. (default False)
```

## Example Run Predict
```
cat alignment_fasta.fa | ./linearalifold
>seq1
CUCACAACGUUUGUGCCUCAGUUACCCGUAGAUGUAGU
>seq2
UCGACACCACU---GCCUCGGUUACCCAUCGGUGCAGU

...((.((((((((((....).....))))))))).)) (-25.30 = 0.00 + -25.31)
```
