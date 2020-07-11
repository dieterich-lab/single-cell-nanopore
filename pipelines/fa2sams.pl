$fa=shift;
$tab=shift;
$adap=shift;
$dir=shift;
$cdna='GCCTGGCTTGTTTGCAAAGGCCCTGGCCAACG';
$clen=length($cdna);
$polyT='T'x20;

open(F,"<$tab");
while(<F>){chomp;
@t=split(/\t/);
$h{$t[0]}=$t[1]
}
close F;

open(F,"<$fa");
while(<F>){chomp;
if(/^>/){@t=split(/\t/,substr($_,1));next}
$_=$adap.$_.$polyT.$cdna;
$l=length($_);
$l2=$l-$clen;
$q='F'x$l;
$h1{$h{$t[0]}}.= "$t[0]..$t[1]\t0\t1\t21\t31\t${l2}S${clen}M\t*\t0\t0\t$_\t$q\n"
}
close F;

foreach(keys %h1){
open(F,">$dir/$_.txt");
print F '@SQ	SN:1	LN:100',"\n";
print F $h1{$_};
close F;
}
