import pandas as pd

# =========================
# 配置参数
# =========================
EXCEL_FILE = '/haplox/users/xuliu/TCR_Project/SampleList/data_information.xlsx'
SHEET_NAME = 'a1'
SERVER_TSV = '/haplox/users/guowei/TCR_Project/LargeCohortList/LN_clonotypefile_meta.tsv'
OUTPUT_FILE = '/haplox/users/xuliu/TCR_Project/SampleList/file2022_complete.tsv'
HEADER_ROW = 2


def main():
    try:
        # -------------------------------------------------
        # 1. 读取 Excel
        # -------------------------------------------------
        df_excel = pd.read_excel(EXCEL_FILE, sheet_name=SHEET_NAME, header=HEADER_ROW)

        # 打印列名（调试时非常有用）
        print("Excel 列名：")
        print(df_excel.columns.tolist())

        # ------------------------------
        # 精确识别列名（不再模糊匹配）
        # ------------------------------
        center_col = next(
            (c for c in df_excel.columns if str(c).strip().lower() == 'center'),
            None
        )

        sample_col = next(
            (c for c in df_excel.columns
             if 'sample' in str(c).lower() and 'id' in str(c).lower()),
            None
        )

        # ⚠️ 关键点：只接受“精确的 Type 列”
        type_col = next(
            (c for c in df_excel.columns if str(c).strip().lower() == 'type'),
            None
        )

        if not center_col:
            raise ValueError("Excel 中未找到 Center 列")
        if not sample_col:
            raise ValueError("Excel 中未找到 Sample_ID 列")
        if not type_col:
            raise ValueError("Excel 中未找到【Type】列（注意不是 Sample_Type）")

        print(f"使用列：Center={center_col}, Sample_ID={sample_col}, Type={type_col}")

        # -------------------------------------------------
        # 2. 清洗关键列（这是防翻车核心）
        # -------------------------------------------------
        df_excel[center_col] = (
            df_excel[center_col]
            .astype(str)
            .str.strip()
            .str.upper()
        )

        df_excel[type_col] = (
            df_excel[type_col]
            .astype(str)
            .str.strip()
            .str.upper()
        )

        df_excel[sample_col] = (
            df_excel[sample_col]
            .astype(str)
            .str.strip()
        )

        # -------------------------------------------------
        # 3. 筛选 Center == SZ 且 Type == BENIGN
        # -------------------------------------------------
        filtered_excel = df_excel[
            (df_excel[center_col] == 'SZ') &
            (df_excel[type_col] == 'BENIGN')
        ]

        print(f"Excel 中筛选得到 {len(filtered_excel)} 个 SZ + Benign 样本")

        if filtered_excel.empty:
            print("⚠️ 提示：请确认 Type 列中是否真的存在 'Benign'")
            print("Type 列实际取值示例：")
            print(df_excel[type_col].value_counts(dropna=False))
            return

        sz_benign_samples = set(filtered_excel[sample_col])

        # -------------------------------------------------
        # 4. 读取 TSV 并匹配
        # -------------------------------------------------
        df_tsv = pd.read_csv(SERVER_TSV, sep='\t')

        tsv_sample_col = next(
            (c for c in df_tsv.columns
             if 'sample' in str(c).lower() and 'id' in str(c).lower()),
            None
        )

        file_col = next(
            (c for c in df_tsv.columns
             if 'file' in str(c).lower() or 'path' in str(c).lower()),
            None
        )

        if not tsv_sample_col or not file_col:
            raise ValueError("TSV 中未找到 Sample_ID 或 File_Path 列")

        df_tsv[tsv_sample_col] = (
            df_tsv[tsv_sample_col]
            .astype(str)
            .str.strip()
        )

        matched_df = df_tsv[df_tsv[tsv_sample_col].isin(sz_benign_samples)]

        # -------------------------------------------------
        # 5. 输出结果
        # -------------------------------------------------
        if matched_df.empty:
            print("未在 TSV 中找到匹配样本")
            return

        result = matched_df[[tsv_sample_col, file_col]].copy()
        result.columns = ['Sample_ID', 'File_Path']

        result.to_csv(OUTPUT_FILE, sep='\t', index=False)

        print(f"成功保存 {len(result)} 条记录到 {OUTPUT_FILE}")

    except Exception as e:
        print(f"错误: {e}")


if __name__ == "__main__":
    main()
