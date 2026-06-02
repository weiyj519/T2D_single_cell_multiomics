#!/usr/bin/env python3

"""
Script name: 02_convert_atac_csv_to_parquet.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/01_feature_matrix_preparation/02_convert_atac_csv_to_parquet.py
"""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src" / "utils"))
from config import ensure_dir, load_config, resolve_path

CONFIG = load_config()
RAW_DATA_DIR = resolve_path(CONFIG, "paths.raw_data_dir")
PROCESSED_DATA_DIR = resolve_path(CONFIG, "paths.processed_data_dir")
SEURAT_OBJECT_PATH = resolve_path(CONFIG, "paths.seurat_object_path")
RNA_EXPRESSION_PARQUET = resolve_path(CONFIG, "paths.RNA_EXPRESSION_PARQUET")
ATAC_PEAK_PARQUET = resolve_path(CONFIG, "paths.ATAC_PEAK_PARQUET")
RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.result_dir"))
FIGURE_DIR = ensure_dir(resolve_path(CONFIG, "paths.figure_dir"))
SUPPLEMENTARY_TABLE_DIR = ensure_dir(resolve_path(CONFIG, "paths.supplementary_table_dir"))
MODEL_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.model_result_dir"))
RNA_MODEL_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.rna_model_result_dir"))
ATAC_MODEL_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.atac_model_result_dir"))
STAGE_GENE_DIR = ensure_dir(resolve_path(CONFIG, "paths.stage_gene_dir"))
STAGE_PEAK_DIR = ensure_dir(resolve_path(CONFIG, "paths.stage_peak_dir"))
DOWNSTREAM_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.downstream_result_dir"))
ATAC_REGULATORY_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.atac_regulatory_result_dir"))
COMMUNICATION_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.communication_result_dir"))
PSEUDOTIME_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.pseudotime_result_dir"))

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

def csv_to_parquet(csv_path, parquet_path=None, compression='snappy', chunk_size=50000):
    """
     CSV  Parquet 
    
    Parameters:
    -----------
    csv_path : str
         CSV 
    parquet_path : str, optional
         Parquet ， None 
    compression : str
        ， 'snappy', 'gzip', 'brotli', 'zstd'
    chunk_size : int
        ，
    """
    csv_path = Path(csv_path)
    
    if parquet_path is None:
        parquet_path = csv_path.with_suffix('.parquet')
    else:
        parquet_path = Path(parquet_path)
    
    print(f"Converting CSV to Parquet...")
    print(f"Input:  {csv_path}")
    print(f"Output: {parquet_path}")
    
    # 
    print("\n1. Reading sample to infer schema...")
    sample_df = pd.read_csv(csv_path, nrows=1000)
    print(f"   Columns: {len(sample_df.columns)}")
    print(f"   Sample rows: {len(sample_df)}")
    
    # 
    print("\n2. Column types:")
    for col in sample_df.columns[:10]:  # 10
        print(f"   {col}: {sample_df[col].dtype}")
    if len(sample_df.columns) > 10:
        print(f"   ... and {len(sample_df.columns) - 10} more columns")
    
    # 
    print("\n3. Reading full CSV file...")
    df = pd.read_csv(csv_path)
    print(f"   Total rows: {len(df)}")
    print(f"   Total columns: {len(df.columns)}")
    print(f"   Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.2f} MB")
    
    #  string （）
    print("\n4. Converting data types...")
    for col in df.columns:
        if df[col].dtype == 'object':
            df[col] = df[col].astype(str)
    
    #  Parquet
    print(f"\n5. Writing to Parquet (compression={compression})...")
    df.to_parquet(
        parquet_path,
        engine='pyarrow',
        compression=compression,
        index=False
    )
    
    # 
    parquet_size = parquet_path.stat().st_size / 1024**2
    print(f"\n✓ Conversion completed!")
    print(f"  Parquet file size: {parquet_size:.2f} MB")
    
    # 
    print("\n6. Verifying Parquet file...")
    test_df = pd.read_parquet(parquet_path, columns=['celltype'])
    print(f"   Successfully read 'celltype' column: {len(test_df)} rows")
    print(f"   Unique cell types: {test_df['celltype'].nunique()}")
    print(f"   Cell types: {sorted(test_df['celltype'].unique())}")
    
    return parquet_path


if __name__ == "__main__":
    csv_path = str(PROCESSED_DATA_DIR / "atac_peaks_autosomes.csv")
    parquet_path = ATAC_PEAK_PARQUET
    
    try:
        output_path = csv_to_parquet(
            csv_path=csv_path,
            parquet_path=parquet_path,
            compression='zstd'  #  'zstd' 
        )
        print(f"\n✓ Successfully converted to: {output_path}")
        
    except Exception as e:
        print(f"\n✗ Error during conversion: {e}")
        import traceback
        traceback.print_exc()