import pandas as pd
import os
import shutil

# 配置参数
MATCHED_FILE = '/haplox/users/xuliu/TCR_Project/SampleList/file2022_complete.tsv'
ROOT_PATH = '/haplox/users/xuliu/TCR_Project/Results/2022_TCRcell'
SUBDIR_NAME = 'Map_Clone_Analysis'

def process_files():
    """处理文件复制和目录创建"""
    try:
        # 读取数据
        df = pd.read_csv(MATCHED_FILE, sep='\t')
        
        if df.empty:
            print("无数据可处理")
            return
        
        print(f"处理 {len(df)} 条记录")
        
        # 创建目录和复制文件
        success_count = 0
        
        for idx, row in df.iterrows():
            sample_id = str(row.get('Sample_ID', '')).strip()
            source_path = str(row.get('File_Path', '')).strip()
            
            if not sample_id or not source_path:
                continue
            
            # 创建目标路径
            target_dir = os.path.join(ROOT_PATH, sample_id, SUBDIR_NAME)
            os.makedirs(target_dir, exist_ok=True)
            
            # 复制文件
            if os.path.exists(source_path):
                filename = os.path.basename(source_path)
                target_path = os.path.join(target_dir, filename)
                
                try:
                    shutil.copy2(source_path, target_path)
                    print(f"✓ {sample_id}: 复制成功")
                    success_count += 1
                except Exception as e:
                    print(f"✗ {sample_id}: 复制失败 - {e}")
            else:
                print(f"✗ {sample_id}: 源文件不存在 - {source_path}")
        
        print(f"\n完成! 成功处理 {success_count}/{len(df)} 个样本")
        
    except Exception as e:
        print(f"错误: {e}")

if __name__ == "__main__":
    process_files()