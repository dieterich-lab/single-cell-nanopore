$fa=shift;
$tab=shift;

open(F,"<$tab");
while(<F>){chomp;
@t=split(/\t/);
$h{$t[0]}=$t[1]
}
close F;

open(F,"<$fa");
while(<F>){
if(/^>/){chomp;@t=split(/\t/,substr($_,1));next}
$t[0]=~s/_\w+$//;
print "$t[0]\t$t[1]\t$h{$t[0]}\t$_" if defined $h{$t[0]}
}
close F;
