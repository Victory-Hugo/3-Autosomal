
# 查找所有Donor_v_Target.chunkcounts.out文件并合并到Chunkcounts.txt，首行只保留一次
outdir="/home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/OUTPUT_bak/东亚受到其他所有"
outfile1="$outdir/Chunkcounts.txt"
first=1
> "$outfile1"
find "$outdir" -type f -name 'Donor_v_Target.chunkcounts.out' | sort | while read file; do
    if [ $first -eq 1 ]; then
        head -n 1 "$file" >> "$outfile1"
        first=0
    fi
    tail -n +2 "$file" >> "$outfile1"
done


outfile2="$outdir/Chunklengths.txt"
first=1
> "$outfile2"
find "$outdir" -type f -name 'Donor_v_Target.chunklengths.out' | sort | while read file; do
    if [ $first -eq 1 ]; then
        head -n 1 "$file" >> "$outfile2"
        first=0
    fi
    tail -n +2 "$file" >> "$outfile2"
done
