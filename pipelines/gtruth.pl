$fa=shift;
$alen=shift;
$blen=shift;
$len=shift;

while(<>){
$h{"$2$3"}="$1\t$2\t$3" if /GN:Z:(\S+).*CB:Z:([ACGT]+).*UB:Z:([ACGT]+)/
}

open(F,"<$fa");
while(<F>){
$b=/^>/;
$s=$1 if /^>(\w+):/;
print $1,"\t",$h{substr($_,$alen,$blen)},"\n" if !$b and defined $h{substr($_,$alen,$blen)}
}
close F
