mkdir('src');
while(<>){
if (/^## (\S+)/){
open(O,">src/$f");
print O $s;
close O;
$f=$1;$b=0;$s=''}
$b++ if /^```/;
$s.=$_ if $b==2;
$b++ if /^```/;
}
