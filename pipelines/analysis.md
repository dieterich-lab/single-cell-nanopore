# gene_names.sh
```
sv possorted_genome_bam.bam|perl -ne 'print "$2\t$1\n" if /GN:Z:(\S+).*CB:Z:([ACGT]+)/' > bc.gene
```
