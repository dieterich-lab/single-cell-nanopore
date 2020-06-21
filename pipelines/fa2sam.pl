$fa=shift;
$tab=shift;
$len=shift;
$clen=shift;
$gene=shift;

open(F,"<$tab");
while(<F>){chomp;
@t=split(/\t/);
$h{$t[0]}="GE:Z:$t[1]\tXX:Z:$t[2]\tUU:Z:$t[3]"
}
close F;

print '@SQ	SN:1	LN:100',"\n";
open(F,"<$fa");
while(<F>){chomp;
$b=/^>/;
$l=length($_);
$l2=$l-$clen;
$q='F'x$l;
$s = substr($_,1) if $b;
$tg = defined $h{$s} ? $h{$s} : "GE:Z:$gene";
print "$s\t0\t1\t21\t31\t${l2}S${clen}M\t*\t0\t0\t$_\t$q\t$tg\tXF:Z:CODING\tUQ:Z:$q\tUS:Z:$_\n" if $l==$len and !$b
}
close F
