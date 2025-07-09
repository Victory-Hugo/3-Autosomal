#!/bin/bash

# 设置错误处理
set -e

# 获取脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# ===== 配置参数（之前放在YAML文件中）=====
PLINK_PATH="plink"
SCRIPTS_DIR="$BASE_DIR/src"
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data"
PLINK2TREEMIX_PATH="plink2treemix.py"
CPU_USAGE=0.8

# ===== 默认值 =====
PLINK_PREFIX="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data"

# 显示帮助信息
show_help() {
    echo "使用方法: $0 [选项] [plink文件前缀]"
    echo "选项:"
    echo "  -p, --plink      PLINK路径 (默认: plink)"
    echo "  -o, --output     输出目录 (默认: $OUTPUT_DIR)"
    echo "  -c, --cpu        CPU使用率 (默认: $CPU_USAGE)"
    echo "  -h, --help       显示此帮助信息"
    echo "示例:"
    echo "  $0 -o /path/to/output /path/to/plink_file"
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--plink)
            PLINK_PATH="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -c|--cpu)
            CPU_USAGE="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            if [[ -z "$PLINK_PREFIX" || "$PLINK_PREFIX" == "/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data" ]]; then
                PLINK_PREFIX="$1"
            fi
            shift
            ;;
    esac
done

# 检查PLINK文件是否存在
if [ ! -f "${PLINK_PREFIX}.bed" ]; then
    echo "错误: PLINK文件 ${PLINK_PREFIX}.bed 不存在"
    exit 1
fi

# 确保输出目录存在
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "无法进入输出目录: $OUTPUT_DIR"; exit 1; }

# 确保Python环境正确
if ! command -v python3 &> /dev/null; then
    echo "错误: 需要Python 3"
    exit 1
fi

# 检查必要的Python包
echo "检查Python依赖..."
python3 -c "import pandas" 2>/dev/null || { echo "错误: 需要pandas包. 运行: pip install pandas"; exit 1; }
python3 -c "import tqdm" 2>/dev/null || { echo "错误: 需要tqdm包. 运行: pip install tqdm"; exit 1; }

# 创建临时Python脚本
TEMP_PY_SCRIPT=$(mktemp)
cat > "$TEMP_PY_SCRIPT" << 'PYTHON_SCRIPT'
#!/usr/bin/env python3

import os
import sys
import subprocess
import multiprocessing as mp
from pathlib import Path
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Tuple, Optional
import pandas as pd
from tqdm import tqdm
import logging

class FSTPipeline:
    def __init__(self, plink_prefix, scripts_dir, output_dir, plink_path, plink2treemix_path, cpu_usage):
        """初始化FST计算流程"""
        self.plink_prefix = Path(plink_prefix)
        self.scripts_dir = Path(scripts_dir)
        self.output_dir = Path(output_dir)
        self.plink_path = plink_path
        self.plink2treemix_path = plink2treemix_path
        self.temp_dir = Path(f"temp_fst_{os.getpid()}")
        self.max_workers = max(1, int(mp.cpu_count() * float(cpu_usage)))
        self.logger = self._setup_logger()
        
    def _setup_logger(self) -> logging.Logger:
        """设置日志"""
        logger = logging.getLogger('FSTPipeline')
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        return logger

    def setup(self):
        """设置工作环境"""
        self.logger.info(f"将使用 {self.max_workers} 个CPU核心进行并行计算")
        os.makedirs(self.temp_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        
    def generate_clusters(self):
        """生成群组文件"""
        self.logger.info("Step 1: 生成群组文件...")
        cluster_script = self.scripts_dir / "makecluster.sh"
        subprocess.run([
            "bash", str(cluster_script),
            str(self.plink_prefix)
        ], check=True)
        
    def read_populations(self) -> List[str]:
        """读取群组信息"""
        self.logger.info("Step 2: 读取群组信息...")
        with open("group.csv") as f:
            return [line.strip() for line in f if line.strip()]
            
    def generate_population_pairs(self, populations: List[str]) -> List[Tuple[str, str]]:
        """生成群组对"""
        self.logger.info("Step 3: 生成群组对...")
        return list(combinations(sorted(populations), 2))
        
    def calculate_fst(self, pair_info: Tuple[int, Tuple[str, str]]) -> Tuple[str, Optional[float]]:
        """计算单个群组对的FST值"""
        idx, (pop1, pop2) = pair_info
        pair_str = f"{pop1} {pop2}"
        temp_prefix = self.temp_dir / f"temp_{idx}"
        log_file = temp_prefix.with_suffix('.log')
        
        try:
            subprocess.run([
                self.plink_path,
                "--bfile", str(self.plink_prefix),
                "--within", f"{str(self.plink_prefix)}.csv",
                "--keep-cluster-names", pair_str,
                "--fst",
                "--out", str(temp_prefix)
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            
            with open(log_file) as f:
                for line in f:
                    if "Weighted" in line:
                        return pair_str, float(line.replace("Weighted Fst estimate: ", "").strip())
        except Exception as e:
            self.logger.warning(f"计算{pair_str}的FST值时出错: {e}")
        
        return pair_str, None
        
    def parallel_fst_calculation(self, population_pairs: List[Tuple[str, str]]) -> pd.DataFrame:
        """并行计算所有群组对的FST值"""
        self.logger.info("Step 4: 并行计算群组间FST值...")
        results = []
        
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            futures = [
                executor.submit(self.calculate_fst, (i, pair))
                for i, pair in enumerate(population_pairs)
            ]
            
            for future in tqdm(as_completed(futures), total=len(futures), desc="计算FST"):
                pair_str, fst_value = future.result()
                if fst_value is not None:
                    pop1, pop2 = pair_str.split()
                    results.append({'Pop1': pop1, 'Pop2': pop2, 'Fst': fst_value})
                
        return pd.DataFrame(results)
        
    def save_results(self, df: pd.DataFrame):
        """保存FST计算结果"""
        self.logger.info("Step 5: 保存FST结果...")
        output_file = self.output_dir / "Fst_mat.csv"
        df.to_csv(output_file, index=False)
        
    def prepare_treemix(self):
        """准备TreeMix输入文件"""
        self.logger.info("Step 6: 准备TreeMix输入文件...")
        freq_file = Path("plink.frq.strat")
        treemix_input = Path("TreeMix.gz")
        
        # 计算频率
        subprocess.run([
            self.plink_path,
            "--bfile", str(self.plink_prefix),
            "--within", f"{str(self.plink_prefix)}.csv",
            "--freq"
        ], check=True)
        
        # 压缩频率文件
        if freq_file.with_suffix('.gz').exists():
            os.remove(str(freq_file.with_suffix('.gz')))
        subprocess.run(["gzip", str(freq_file)], check=True)
        
        # 运行plink2treemix.py
        plink2treemix_full_path = self.scripts_dir / self.plink2treemix_path
        subprocess.run([
            str(plink2treemix_full_path),
            str(freq_file.with_suffix('.gz')),
            str(treemix_input)
        ], check=True)
        
    def cleanup(self):
        """清理临时文件"""
        self.logger.info("Step 7: 清理临时文件...")
        if self.temp_dir.exists():
            subprocess.run(["rm", "-rf", str(self.temp_dir)])
        for pattern in ["*.P", "*.Q"]:
            for f in Path().glob(pattern):
                try:
                    f.unlink()
                except:
                    pass
                
    def run(self):
        """运行完整的分析流程"""
        try:
            self.setup()
            self.generate_clusters()
            populations = self.read_populations()
            population_pairs = self.generate_population_pairs(populations)
            results_df = self.parallel_fst_calculation(population_pairs)
            self.save_results(results_df)
            self.prepare_treemix()
            self.cleanup()
            self.logger.info("分析完成！")
        except Exception as e:
            self.logger.error(f"错误: {str(e)}")
            import traceback
            traceback.print_exc()
            self.cleanup()
            sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(f"Usage: {sys.argv[0]} <plink_prefix> <scripts_dir> <output_dir> <plink_path> <plink2treemix_path> <cpu_usage>")
        sys.exit(1)
    
    pipeline = FSTPipeline(
        plink_prefix=sys.argv[1],
        scripts_dir=sys.argv[2],
        output_dir=sys.argv[3],
        plink_path=sys.argv[4],
        plink2treemix_path=sys.argv[5],
        cpu_usage=sys.argv[6]
    )
    pipeline.run()
PYTHON_SCRIPT

echo "开始运行分析..."
cd "$OUTPUT_DIR" || { echo "无法进入输出目录: $OUTPUT_DIR"; exit 1; }

# 运行Python脚本
python3 "$TEMP_PY_SCRIPT" \
    "$PLINK_PREFIX" \
    "$SCRIPTS_DIR" \
    "$OUTPUT_DIR" \
    "$PLINK_PATH" \
    "$PLINK2TREEMIX_PATH" \
    "$CPU_USAGE"

# 清理临时Python脚本
rm -f "$TEMP_PY_SCRIPT"

echo "分析完成！"