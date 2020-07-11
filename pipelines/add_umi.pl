$tab=shift;
$umi=shift;

open(F,"<$umi");
while(<F>){chomp;
@t=split(/\t/);
$h{$t[0]}=1
}
close F;

open(F,"<$tab");
$_=<F>;print;
while(<F>){chomp;
@t=split(/\t/);
$t[$#t]=1 if defined $h{"$t[0]..$t[1]"};
print join "\t", @t, "\n"; 
}
close F;

