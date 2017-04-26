![logo](http://www.bi.cs.titech.ac.jp/ghostz-gpu/ghostzgpu-logo.png)
======

GHOSTZ-GPU is a homology search tool which can detect remote homologues like BLAST and is about 5-7 times more efficient than GHOSTZ by using GHOSTZ. 

GHOSTZ-GPU outputs search results in the format similar to BLAST-tabular format.

Requirements
------------
- gcc => 4.3
- Boost >= 1.55.0
- CUDA >= 6.0

Installation
------------
1. Download the archive of GHOSTZ-GPU from this repository.
2. Extract the archive and cd into the extracted directory.
3. Run make command.
4. Copy `ghostz-gpu` binary file to any directory you like.

Commands:

    $ tar xvzf ghostz-gpu.tar.gz
    $ cd ghostz-gpu
    $ make BOOST_PATH=Boost CUDA_TOOLKIT_PATH=CUDA
    $ cp ghostz-gpu /AS/YOU/LIKE/


Boost and CUDA are directories where they are installed, respectively.
    
Usage
-----
GHOSTZ-GPU requires specifically formatted database files for homology search. These files can be generated from FASTA formatted DNA/protein sequence files. 

Users have to prepare a database file in FASTA format and convert it into GHOSTZ-GPU format database files by using GHOSTZ-GPU `db` command at first. GHOSTZ-GPU `db` command requires 2 args (`[-i dbFastaFile]` and `[-o dbName]`). GHOSTZ-GPU `db` command divides a database FASTA file into several database chunks and generates several files (.inf, .ind, .nam, .pos, .seq). All generated files are needed for the search. Users can specify the size of each chunk. Smaller chunk size requires smaller memory, but efficiency of the search will decrease.  
For executing homology search, GHOSTZ-GPU `aln` command is used and that command requires at least 2 args (`[-i qryName]` and `[-d dbName]`).

Example
-------

    $ ghostz-gpu db  -i ./data/db.fasta -o exdb
    $ ghostz-gpu aln -q d -t p -i ./data/queries.fasta -d exdb -o exout

Command and Options
-------------------
`db`: convert a FASTA file to GHOSTZ format database files

    ghostz-gpu db [-i dbFastaFile] [-o dbName] [-C clustering][-l chunkSize]
        [-L clusteringSubsequenceLength]  [-s seedThreshold]
    Options:
    (Required)
      -i STR    Protein sequences in FASTA format for a database
      -o STR    The name of the database
    (Optional)
      -C STR    Clustering, T (enable) or F (disable) [T]
      -l INT    Chunk size of the database (bytes) [1073741824 (=1GB)]
      -L INT    Length of a subsequence for clustering [10]
      -s INT    The seed threshold [39]


`aln`:  Search homologues of queries from database

    ghostz-gpu aln [-i queries] [-o output] [-d database] [-v maxNumAliSub]
      [-b maxNumAliQue] [-h hitsSize] [-l queriesChunkSize] [-q queryType]
      [-t databaseType] [-F filter] [-a numThreads] [-g numGPUs]
    Options:
    (Required)
      -i STR    Sequences in FASTA format
      -o STR    Output file
      -d STR    database name (must be formatted)
    (Optional)
      -v INT    Maximum number of alignments for each subject [1]
      -b INT    Maximum number of the output for a query [10]
      -l INT    Chunk size of the queries (bytes) [134217728 (=128MB)]
      -q STR    Query sequence type, p (protein) or d (dna) [p]
      -t STR    Database sequence type, p (protein) or d (dna) [p]
      -F STR    Filter query sequence, T (enable) or F (disable) [T] 
      -a INT    The number of threads [1]
      -g INT    The number of GPUs [the number of available GPUs]

Search results
--------------
GHOSTZ-GPU outputs the tab-deliminated file as search results.

Example)

```
query0  subject0        100     25      0       0       1       75      1       25      2.75456e-15     60.4622
query0  subject6        100     10      0       0       46      75      16      25      2.58417e-05     27.335
query1  subject0        100     24      0       0       2       73      1       24      1.36707e-14     58.151
query1  subject6        100     9       0       0       47      73      16      24      0.000128251     25.0238
query3  subject6        100     14      0       0       34      75      12      25      3.60591e-07     33.4982
query3  subject0        100     10      0       0       46      75      16      25      2.58417e-05     27.335
query4  subject6        100     14      0       0       42      1       12      25      3.60591e-07     33.4982
query4  subject0        100     10      0       0       30      1       16      25      2.58417e-05     27.335

```

Each column shows;

1. Name of a query sequence
2. Name of a homologue sequence (subject)
3. Sequence Identity
4. Alignment length
5. The number of mismatches in the alignment
6. The number of gap openingsin the alignemt
7. Start position of the query in the alignment
8. End position of the query in the alignemnt
9. Start position of the subject in the alignment
10. End position of the subject in the alignment
11. *E*-value
12. Normalized score

References
----------
Shuji Suzuki, Masanori Kakuta, Takashi Ishida, Yutaka Akiyama. Acceleration of sequence homology searches by means of graphics processing units (GPUs) and database subsequence clustering. (submitted)


Copyright © 2015 Akiyama_Laboratory, Tokyo Institute of Technology, All Rights Reserved.  

