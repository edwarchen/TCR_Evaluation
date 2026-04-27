import csv
from pathlib import Path

# 原始 primer 文件
input_file = "data/primer_set/tcr引物序列.xlsx"
v_output = "data/TRBV_primers.fasta"
j_output = "data/TRBJ_primers.fasta"


def iter_primer_rows(file_path: Path):
    """Yield (primer_id, sequence) from csv/xlsx first two columns."""
    suffix = file_path.suffix.lower()

    if suffix in {".xlsx", ".xls"}:
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "读取 Excel 需要 pandas 和 openpyxl，请先安装：pip install pandas openpyxl"
            ) from exc

        df = pd.read_excel(file_path, engine="openpyxl")
        for _, row in df.iloc[:, :2].iterrows():
            primer_id = str(row.iloc[0]).strip() if len(row) > 0 else ""
            sequence = str(row.iloc[1]).strip() if len(row) > 1 else ""
            yield primer_id, sequence
        return

    if suffix == ".csv":
        with open(file_path, "r", encoding="utf-8-sig", newline="") as f_in:
            reader = csv.reader(f_in)
            for row in reader:
                if not row or len(row) < 2:
                    continue
                yield row[0].strip(), row[1].strip()
        return

    raise ValueError(f"不支持的文件类型: {file_path.suffix}")


def main():
    src = Path(input_file)
    if not src.exists():
        print(f"错误: 找不到文件 '{input_file}'，请确认文件路径是否正确。")
        return

    v_count = 0
    j_count = 0

    with open(v_output, "w", encoding="utf-8") as f_v, open(
        j_output, "w", encoding="utf-8"
    ) as f_j:
        for primer_id, sequence in iter_primer_rows(src):
            if not primer_id or not sequence:
                continue
            if primer_id.upper() == "ID" or sequence.upper() == "SEQUENCE":
                continue

            # 将斜杠和空格替换为下划线，避免下游工具解析异常
            clean_id = primer_id.replace("/", "_").replace(" ", "_")

            if "V" in clean_id.upper():
                f_v.write(f">{clean_id}_fw\n{sequence}\n")
                v_count += 1
            elif "J" in clean_id.upper():
                f_j.write(f">{clean_id}_rev\n{sequence}\n")
                j_count += 1
            else:
                print(f"警告: 无法识别引物方向，已跳过 -> {primer_id}")

    print("处理完成")
    print(f"成功提取了 {v_count} 条 TRBV 正向引物，已保存至: {v_output}")
    print(f"成功提取了 {j_count} 条 TRBJ 反向引物，已保存至: {j_output}")


if __name__ == "__main__":
    main()
