#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;
use List::Util qw(max);

# 获取命令行参数
my $fileName = shift @ARGV or die "Usage: $0 <plink_file_prefix>\n";

# 定义变量
my @populations;
my @population_pairs;
my @fst_results;

# 设置并行进程数（根据CPU核心数自动设置）
my $cpu_cores = `nproc`;
chomp($cpu_cores);
my $max_processes = max(1, int($cpu_cores * 0.8));  # 使用80%的CPU核心
print "将使用 $max_processes 个CPU核心进行并行计算\n";

# 创建临时目录用于并行处理
my $temp_dir = "temp_fst_$$";
system("mkdir -p $temp_dir");

# 第一步：生成群组文件
print "Step 1: 生成群组文件...\n";
system("bash /mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/1-treemix/src/makecluster.sh $fileName");

# 第二步：读取群组信息
print "Step 2: 读取群组信息...\n";
open(my $group_file, "<", "group.csv") or die "无法打开group.csv: $!";
while (my $line = <$group_file>) {
    chomp($line);
    push(@populations, $line);
}
close($group_file);

# 第三步：生成群组对
print "Step 3: 生成群组对...\n";
my %seen_pairs;
foreach my $pop1 (@populations) {
    foreach my $pop2 (@populations) {
        next if $pop1 eq $pop2;
        my $pair = join(" ", sort($pop1, $pop2));
        next if $seen_pairs{$pair}++;
        push(@population_pairs, "$pop1 $pop2");
    }
}

# 第四步：并行计算FST
print "Step 4: 并行计算群组间FST值...\n";
my $pm = Parallel::ForkManager->new($max_processes);

# 设置数据收集回调
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
    if (defined($data_ref)) {
        push(@fst_results, $$data_ref);
    }
});

foreach my $pair_idx (0 .. $#population_pairs) {
    my $pair = $population_pairs[$pair_idx];
    
    # 启动子进程
    $pm->start and next;
    
    # 在子进程中计算FST
    my $temp_log = "$temp_dir/plink_$pair_idx.log";
    system("plink --bfile $fileName --within $fileName.csv --keep-cluster-names '$pair' --fst --out $temp_dir/temp_$pair_idx 2> $temp_log");
    
    # 提取FST值
    my $result;
    open(my $fst_file, "<", "$temp_dir/temp_$pair_idx.log") or die "无法打开 $temp_dir/temp_$pair_idx.log: $!";
    while (my $line = <$fst_file>) {
        if ($line =~ /Weighted/) {
            $result = "$pair,$line";
            last;
        }
    }
    close($fst_file);
    
    $pm->finish(0, \$result);  # 发送结果回主进程
}

$pm->wait_all_children;  # 等待所有子进程完成

# 第五步：保存FST结果
print "Step 5: 保存FST结果...\n";
open(my $out_file, ">", "fst.csv") or die "无法创建fst.csv: $!";
print $out_file join("", @fst_results);
close($out_file);

# 清理和格式化FST结果
system("sed -i.bu 's/Weighted Fst estimate: //g' fst.csv");
system("sed -i.bu 's/,/ /g' fst.csv");
system("echo 'Pop1 Pop2 Fst' > Fst.txt");
system("cat Fst.txt fst.csv > Fst_mat.csv");

# 第六步：准备TreeMix输入文件
print "Step 6: 准备TreeMix输入文件...\n";
system("plink --bfile $fileName --within $fileName.csv --freq");
system("rm -f plink.frq.strat.gz") if (-e "plink.frq.strat.gz");
system("gzip plink.frq.strat");
system("/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/1-treemix/src/plink2treemix.py plink.frq.strat.gz TreeMix.gz");

# 清理临时文件
print "Step 7: 清理临时文件...\n";
system("rm -rf $temp_dir");
system("rm -f *.P *.Q");

print "分析完成！\n";

