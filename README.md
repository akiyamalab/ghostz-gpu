GHOSTZ
======

GHOSTZ is a homology search tool which can detect remote homologues like BLAST and is about 200 times more efficient than BLAST by using database subsequence clustering. GHOSTZ outputs search results in the format similar to BLAST-tabular format.

Requirements
------------
gcc >4.3

Installation
------------

Download the archive of GHOSTZ from the above link.
Extract the archive and cd into the extracted directory.
Run make command.
Copy 'ghostz' file to any directory you like.
      $ tar xvzf ghostz.tar.gz
      $ cd ghostz
      $ make
      $ cp ghostz /AS/YOU/LIKE/
    
Usage
-----

GHOSTZ requires specifically formatted database files for homology search. These files can be generated from FASTA formatted DNA/protein sequence files. 
Users have to prepare a database file in FASTA format and convert it into GHOSTZ format database files by using GHOSTZ "db" command at first. GHOSTZ"db" command requires 2 args ([-i dbFastaFile] and [-o dbName]). GHOSTZ "db" command divides a database FASTA file into several database chunks and generates several files (.inf, .ind, .nam, .pos, .seq). All generated files are needed for the search. Users can specify the size of each chunk. Smaller chunk size requires smaller memory, but efficiency of the search will decrease. 
For executing homology search, GHOSTZ "aln" command is used and that command requires at least 2 args([-i qryName] and [-d dbName]).

Example
-------
$ ghostz db  -i ./db.fasta -o exdb

$ ghostz aln -i exqry -d exdb -o exout

Command and Options
-------------------
db: convert a FASTA file to GHOSTX format database files

  ghostz db [-i dbFastaFile] [-o dbName] [-C clustering][-l chunkSize]
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


aln:  Search homologues of queries from database

  ghostz aln [-i queries] [-o output] [-d database] [-v maxNumAliSub]
             [-b maxNumAliQue] [-h hitsSize] [-l queriesChunkSize] [-q queryType]
             [-t databaseType] [-F filter] [-a numThreads]

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
Search results
--------------
GHOSTZ outputs the tab-deliminated file as search results.

Example)
  hsa:124045...   hsa:124045...   100       139     0       0       1       139     1       139     2.04391e-76     283.878
  hsa:124045...   ptr:454320...   99.2126        127     1       0       13      139     14      140     5.96068e-68     255.758
  hsa:124045...   mcc:714360...   88.9764        127     14      0       13      139     14      140     5.05773e-59     226.098
  hsa:124045...   rno:292078...   58.6777        121     46      2       13      133     14      130     1.38697e-32     138.272
  hsa:124045...   mmu:320869...   55.9055        127     50      3       13      139     12      132     1.17414e-31     135.191
  hsa:124045...   pon:100434...   96.4912        57      2       0       13      69      14      70      3.65839e-25     113.62
  hsa:124045...   bta:100335...   44.9275        138     71      3       2       139     25      157     4.04482e-24     110.153
  hsa:124045...   aml:100464...   26.6667        75      46      2       13      81      1183    1254    0.820692        32.7278
  hsa:124045...   bfo:BRAFLD...   56    25      10      1       108     131     581     605     0.820692        32.7278
  hsa:124045...   tgu:100227...   26.1682        107     69      3       25      130     150     247     1.82831 31.5722

Each column shows;
1.  Name of a query sequence
2.  Name of a homologue sequence (subject)
3.  Sequence Identity
4.  Alignment length
5.  The number of mismatches in the alignment
6.  The number of gap openingsin the alignemt
7.  Start position of the query in the alignment
8.  End position of the query in the alignemnt
9.  Start position of the subject in the alignment
10. End position of the subject in the alignment
11. E-value
12. Normalized score

Copyright © 2014 Akiyama_Laboratory , Tokyo Institute of Technology , All Rights Reserved.  

