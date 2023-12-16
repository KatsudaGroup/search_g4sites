# search_g4sites

## Preparation

`wget` command is essential to download files.

```
$ pip install -r requirements.txt
$ sh download.sh	# Download the fasta files of the human mRNAs.
```

## Analyze fasta files.

```
$ python search_g4sites.py data/human.{1..12}.rna.fna
```

if you want to change the length of the loop length limit, `-n (length)` is available.

```
$ python search_g4sites.py -n 8 data/human.{1..12}.rna.fna
```
