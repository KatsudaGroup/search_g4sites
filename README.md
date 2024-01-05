# search_g4sites

This code is used for the estimation of the number of genes that Staple oligomer approach is potentially applicable.

## Preparation

First of all, we have to prepare the dependencies and download the datasets from the NCBI.
To download the datasets from NCBI, `wget` command is essential.

```
$ python -m venv env; source env/bin/activate	# Optional; create the separated environment.
$ pip install -r requirements.txt				# Install the dependencies
$ sh download.sh	# Download the fasta files of the human mRNAs.
```

## Analyze Fasta files

Analysis is divided into two steps as follows:

1. Count the number of motifs of each mRNA, by using regular expression.
2. Summarize the result.

### Count motifs

```
$ python search_g4sites.py data/human.{1..12}.rna.fna 
```

if you want to change the length of the loop length limit, `-n (length)` is available.

```
$ python search_g4sites.py -n 8 data/human.{1..12}.rna.fna
```

The code prints the number of the motifs in each mRNA. 
In order to summarize, the outputs must be printed in a text file by redirection or `tee` command as follows.

```
$ python search_g4sites.py data/human.{1..12}.rna.fna | tee search_result.tsv
``` 

### Summarize


```
$ python analyze_result.py search_result.tsv
```