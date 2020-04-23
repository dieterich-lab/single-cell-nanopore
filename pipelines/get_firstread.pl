$L=100;
while(<>){
	chomp;
	if (/^>/){
		$id=$_
	}else{
		@t=split(/_/,$id);
		$s=$_;
		$d=$t[5];
		if ($t[4] eq 'R'){
			$s=reverse $s;
			$s=~tr/ATGCatgc/TACGtacg/;
			$d=$t[7]
		}
		$s=substr($s,$d,$t[6]);
		print $id,"\n",$t[6]>$L?substr($s,$L-$t[1]%$L,68):$s,"\n"
	}
}
