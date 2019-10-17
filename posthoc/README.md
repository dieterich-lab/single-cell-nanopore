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
java -cp bk.jar:JACUSA2-ONT.jar qw.Predict -m model.txt -i test_leftTail.log -o output.txt
```
where "output.txt" contains the filtered barcodes. 
Simulate 5000 random reads from the model
```sh
java -cp bk.jar:JACUSA2-ONT.jar qw.Simulate -m model.txt -i cb.fa -o sim.sam
```

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
The output contains three columns: read id, assigned barcode, bayes factor
```
beccdc96-900a-475e-a1e9-6c6c29e5fe68    TTGCTGCAGGCGTCCT        2.828E9
d67416fb-5580-47f4-9113-ef3a3b3a3010    CCCTGATCACCTTCCA        759.535
2cbcce34-d43f-4182-ae30-f8c720fd8a2f    CCACTTGTCGTCAACA        5.36E7
ed454751-bde8-4817-83b2-35a2ffffe8b7    CCGGTGACACTTACAG        387.59
4d204b79-5f84-4e42-b358-2aa9de6e9252    GAGGGTATCGGCTGTG        5.15E7
8b9d6c7c-4566-4ae9-b5d8-b4069a073d5a    TTCCGTGCACATGACT        2.824E9
b25fc275-9eac-419b-b43e-a55a37cb9920    GAACACTCAACAAAGT        2.8E9
```
