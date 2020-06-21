$fa=shift;
$alen=shift;
$blen=shift;
$ulen=shift;
$gene=shift;
$bulen=$blen+$ulen;
$barcode='N'x$blen;
$umi='N'x$ulen;

while(<>){
$h{"$2$3"}="$1\t$2\t$3" if /GN:Z:(\S+).*CB:Z:([ACGT]+).*UB:Z:([ACGT]+)/
}

open(F,"<$fa");
while(<F>){
$b=/^>/; 
$s=$1 if /^>(\w+):/;
$d = defined $h{substr($_,$alen,$bulen)} ? $h{substr($_,$alen,$bulen)} : "$gene\t$barcode\t$umi";
print "$1\t$d\n" unless $b
}
close F
~      
