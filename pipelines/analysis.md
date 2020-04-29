## gene_names.sh
```
sv possorted_genome_bam.bam|perl -ne 'print "$2\t$1\n" if /GN:Z:(\S+).*CB:Z:([ACGT]+)/' > bc.gene
awk -v OFS='\t' '{print $2,$1}' bc.gene > gene.bc
```
## count_umi.sh
```
perl -F"\t" -ane '$h{$F[0]}{$F[1]}++; END { print "$_\t".(keys %{$h{$_}})."\n" for keys %h }' fc1.barumi > fc1.barumi.cnt
```
