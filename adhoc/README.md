# Ad hoc analysis of Flexbar barcode alignment
## Command line wrapper
Parameter estimation 
```sh
f="$1"
java -cp bk.jar:JACUSA_v2.0.0-RC11.jar:commons-cli-1.4.jar:commons-math3-3.6.1.jar:htsjdk-2.12.0-SNAPSHOT.jar qw.Train call-1 -R GRCh38_90.fa -p 2 -c 1 -m 0 -q 2 -P UNSTRANDED -r model.txt $f
```
Filter barcode alignment
```sh
f="$1"
java -cp bk.jar:JACUSA_v2.0.0-RC11.jar:commons-cli-1.4.jar:commons-math3-3.6.1.jar:htsjdk-2.12.0-SNAPSHOT.jar qw.Predict -m model.txt -i $f -o ${f%.*}.barcodelist
```
where *.barcodelist contains the filtered barcodes. 

## Model file specification
A model file contains six blocks as the following example:
```r
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
