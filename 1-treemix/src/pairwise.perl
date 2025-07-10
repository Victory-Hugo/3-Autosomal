#!/usr/bin/perl

($fileName, $treemix_mig_num)=@ARGV;

#!这个文件是通过makecluster.sh生成的。无需处理。
open($file,"group.csv");
@ary;
@ary2;
@ary3;
@ary4;

while ($line=<$file>){
chomp($line);
push(@ary,$line);
push(@ary2,$line);
}#读取参数结束。


foreach $v(@ary)
{
	foreach $v2(@ary2)
	{
	if ($v2 ne $v)
		{
			$invert="$v2 $v";
			unless (grep(/$invert/i,@ary3))
			{
				push(@ary3,"$v $v2");
			}
		}
	}
}

close($file);
# 生成了一个数组，包含了所有的两两组合。
foreach $v(@ary3){
system("rm plink.log");
system( "plink --bfile $fileName --within $fileName.csv --keep-cluster-names $v --fst --threads 16"); 

open($fst,"plink.log");

while($line=<$fst>){

	if($line=~"Weighted"){

		push(@ary4,$v.",".$line);
	}
}

open($fst2,">","fst.csv");
foreach $v(@ary4){

	print $fst2 $v;
}

system("sed -i.bu 's/Weighted Fst estimate: //g' fst.csv");


system("sed -i.bu 's/,/ /g' fst.csv");

system("echo 'Pop1 Pop2 Fst' > Fst.txt ");
system("cat Fst.txt fst.csv > Fst_mat.csv");

}
# for treemix

system("plink --bfile $fileName --within $fileName.csv --freq")

