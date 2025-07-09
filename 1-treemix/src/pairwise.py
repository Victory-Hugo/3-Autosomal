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
import yaml
import logging
from dataclasses import dataclass

@dataclass
class Config:
    """配置数据类"""
    scripts_dir: Path
    output_dir: Path
    plink_path: str
    plink2treemix_path: str
    cpu_usage: float

class FSTPipeline:
    def __init__(self, plink_prefix: str, config: Config):
        """
        初始化FST计算流程
        
        Args:
            plink_prefix: PLINK文件前缀路径
            config: 配置对象
        """
        self.plink_prefix = Path(plink_prefix)
        self.config = config
        self.temp_dir = Path(f"temp_fst_{os.getpid()}")
        self.max_workers = max(1, int(mp.cpu_count() * config.cpu_usage))
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
        os.makedirs(self.config.output_dir, exist_ok=True)
        
    def generate_clusters(self):
        """生成群组文件"""
        self.logger.info("Step 1: 生成群组文件...")
        cluster_script = self.config.scripts_dir / "makecluster.sh"
        subprocess.run([
            "bash", str(cluster_script),
            str(self.plink_prefix)
        ], check=True)
        
    def read_populations(self) -> List[str]:
        """读取群组信息"""
        self.logger.info("Step 2: 读取群组信息...")
        group_file = self.config.output_dir / "group.csv"
        with open(group_file) as f:
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
                self.config.plink_path,
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
        output_file = self.config.output_dir / "Fst_mat.csv"
        df.to_csv(output_file, index=False)
        
    def prepare_treemix(self):
        """准备TreeMix输入文件"""
        self.logger.info("Step 6: 准备TreeMix输入文件...")
        freq_file = self.config.output_dir / "plink.frq.strat"
        treemix_input = self.config.output_dir / "TreeMix.gz"
        
        # 计算频率
        subprocess.run([
            self.config.plink_path,
            "--bfile", str(self.plink_prefix),
            "--within", f"{str(self.plink_prefix)}.csv",
            "--freq",
            "--out", str(freq_file.with_suffix(''))
        ], check=True)
        
        # 压缩频率文件
        if freq_file.with_suffix('.gz').exists():
            freq_file.with_suffix('.gz').unlink()
        subprocess.run(["gzip", str(freq_file)], check=True)
        
        # 运行plink2treemix.py
        subprocess.run([
            str(self.config.scripts_dir / self.config.plink2treemix_path),
            str(freq_file.with_suffix('.gz')),
            str(treemix_input)
        ], check=True)
        
    def cleanup(self):
        """清理临时文件"""
        self.logger.info("Step 7: 清理临时文件...")
        if self.temp_dir.exists():
            subprocess.run(["rm", "-rf", str(self.temp_dir)])
        for pattern in ["*.P", "*.Q"]:
            for f in Path(self.config.output_dir).glob(pattern):
                f.unlink()
                
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
            self.cleanup()
            sys.exit(1)

def load_config(config_path: str) -> Config:
    """加载配置文件"""
    with open(config_path) as f:
        config_data = yaml.safe_load(f)
    
    return Config(
        scripts_dir=Path(config_data['paths']['scripts_dir']),
        output_dir=Path(config_data['paths']['output_dir']),
        plink_path=config_data['tools']['plink'],
        plink2treemix_path=config_data['tools']['plink2treemix'],
        cpu_usage=config_data['compute']['cpu_usage']
    )

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 pairwise.py <config_file> <plink_file_prefix>")
        sys.exit(1)
        
    config = load_config(sys.argv[1])
    pipeline = FSTPipeline(sys.argv[2], config)
    pipeline.run()

if __name__ == "__main__":
    main()
