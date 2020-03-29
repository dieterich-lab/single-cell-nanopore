print '@SQ      SN:1    LN:100',"\n";
while(<>){chomp;
$b=/^>/;
$l=length($_);
$q='I'x$l;
print substr($_,1),"\t" if $b;
print "0\t1\t21\t31\t${l}S\t*\t0\t0\t$_\t$q\n" unless $b;
}

