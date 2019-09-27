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
