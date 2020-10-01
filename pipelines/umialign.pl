$illu=shift;
$nano=shift;
$dir=shift;
$dir='.' unless defined $dir;
$n=shift;

open(F,"<$illu");
while(<F>){chomp;
@t=split(/\t/);
$h{"$t[2]\t$t[0]"}.=">$t[1]\n$t[1]\n"
}
close F;

open(F,"<$nano");
while(<F>){
@t=split(/\t/);
$h1{"$t[2]\t$t[1]"}.=">$t[0]\n$t[3]"
}
close F;

foreach $g(keys %h1){
next unless defined $h{$g};
$r='';
open(F,">$dir/s$n.fa")||die "$!";
print F $h1{$g};
close F;
open(F,">$dir/q$n.fa")||die "$!";
print F $h{$g};
close F;
system("/home/qwang/bin/mummer-4.0.0beta2/mummer q$n.fa s$n.fa -maxmatch -b -c -l 7 -F > mm$n");
open(F,"<$dir/mm$n");
while(<F>){
$b=/^>/;
$ss=/Reverse/?1:0 if $b;
$_=substr($_,2);
$_=substr($_,0,index($_,' '));
$r.="$g\t$_\t$s\t$ss\n" unless $b;
$s=$_ if $b
}
close F;
print $r
}
