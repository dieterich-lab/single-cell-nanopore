# Post hoc analysis of Flexbar barcode alignment
## Command line wrapper
Please find the required sample files from the "test" folder under the repository. 
Parameter estimation 
```sh
java -cp bk.jar:JACUSA2-ONT.jar qw.Train -R ref.fa -p 2 -c 1 -m 2 -q 2 -r model.txt nanopore.bam
java -cp bk.jar:JACUSA2-ONT.jar qw.CB -R ref.fa -p 2 -c 1 -m 2 -q 2 -r model.txt illumina.bam
```
Filter barcode alignment
```sh
java -cp bk.jar:JACUSA2-ONT.jar qw.Predict -m model.txt -i alignment.txt -o output.txt
```
where "output.txt" contains the filtered barcodes. 

## Model file specification
A model file contains six blocks as the following example:
```
#Base substitutions
AA      129806227
AC      507845
AG      1786250
AT      603749
...
#Insertion freq
0       3522983
1       102491
2       47026
...
#Deletion freq
0       3483535
1       118794
2       54391
...
#Insertions
TGAGGTA 10
CCAAGAAA        19
GACAACT 11
...
#Deletions
GGAGCGGC        31
GATTGTGA        15
GATTGTGC        6
N       27451
#Barcodes
AAACCCAAGACACACG        27305
AAACCCATCCCTGGTT        45955
AAACGAAAGGCTGGAT        6823
AAACGCTCATCGGCCA        27690
AAACGCTGTGACAACG        27692
AAACGCTGTTATGTGC        36930
AAAGAACGTTCTCCTG        979
...
```

## Output file format
The output contains three columns: read id, assigned barcode, bayes factor (log)
```
b617bd33-d719-4525-b975-c3f7e7213e53    CGGAACCGTGTACATC        7.0648233761436
beccdc96-900a-475e-a1e9-6c6c29e5fe68    TTGCTGCAGGCGTCCT        21.783081000863707
ed454751-bde8-4817-83b2-35a2ffffe8b7    CCGGTGACACTTACAG        11.953714653236243
f345b646-98dc-48ae-aebf-33fd7fbf8c28    CAGATTGCATGAAAGT        7.077305491988227
4d204b79-5f84-4e42-b358-2aa9de6e9252    GAGGGTATCGGCTGTG        17.75756713842614
ad202c2f-830a-49c9-bc96-7837ac5c8430    GGGCTACGTATAGGGC        7.207479839616639
1272c4a1-7d2a-4988-aa71-d0b4a69c3987    CAGAGCCTCTTGTTAC        17.8281203773929
```
