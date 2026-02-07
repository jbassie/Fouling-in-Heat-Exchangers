"""
Dataset Physics Consistency Validation

This script checks existing datasets for physics consistency:
- Computes correlations between fouling resistance and heat transfer metrics
- Flags rows where fouling increases but Q/U also increase (should be negative correlation)
"""
import pandas as pd
import numpy as np
import os


def validate_dataset_physics(csv_path: str) -> dict:
    """
    Validate physics consistency of a PHE dataset.
    
    Args:
        csv_path: Path to the CSV dataset file
        
    Returns:
        dict: Validation results with correlations and flags
    """
    if not os.path.exists(csv_path):
        print(f"Error: File not found: {csv_path}")
        return {}
    
    print(f"Loading dataset: {csv_path}")
    df = pd.read_csv(csv_path)
    
    print(f"Dataset shape: {df.shape}")
    print(f"Columns: {list(df.columns)}\n")
    
    # Check required columns
    required_cols = ['R_f_hot', 'R_f_cold', 'U_overall', 'Q_actual']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        return {}
    
    results = {
        'file': csv_path,
        'n_rows': len(df),
        'correlations': {},
        'inconsistencies': [],
    }
    
    # Calculate correlations
    print("=== Correlation Analysis ===")
    
    # R_f_hot vs U_overall (should be negative)
    corr_hot_U = df['R_f_hot'].corr(df['U_overall'])
    results['correlations']['R_f_hot_vs_U'] = corr_hot_U
    print(f"R_f_hot vs U_overall: {corr_hot_U:.4f} (expected: negative)")
    
    # R_f_hot vs Q_actual (should be negative)
    corr_hot_Q = df['R_f_hot'].corr(df['Q_actual'])
    results['correlations']['R_f_hot_vs_Q'] = corr_hot_Q
    print(f"R_f_hot vs Q_actual: {corr_hot_Q:.4f} (expected: negative)")
    
    # R_f_cold vs U_overall (should be negative)
    corr_cold_U = df['R_f_cold'].corr(df['U_overall'])
    results['correlations']['R_f_cold_vs_U'] = corr_cold_U
    print(f"R_f_cold vs U_overall: {corr_cold_U:.4f} (expected: negative)")
    
    # R_f_cold vs Q_actual (should be negative)
    corr_cold_Q = df['R_f_cold'].corr(df['Q_actual'])
    results['correlations']['R_f_cold_vs_Q'] = corr_cold_Q
    print(f"R_f_cold vs Q_actual: {corr_cold_Q:.4f} (expected: negative)")
    
    # Total fouling vs metrics
    if 'R_f_total' not in df.columns:
        df['R_f_total'] = df['R_f_hot'] + df['R_f_cold']
    
    corr_total_U = df['R_f_total'].corr(df['U_overall'])
    corr_total_Q = df['R_f_total'].corr(df['Q_actual'])
    results['correlations']['R_f_total_vs_U'] = corr_total_U
    results['correlations']['R_f_total_vs_Q'] = corr_total_Q
    print(f"R_f_total vs U_overall: {corr_total_U:.4f} (expected: negative)")
    print(f"R_f_total vs Q_actual: {corr_total_Q:.4f} (expected: negative)")
    
    # Flag inconsistencies
    print("\n=== Inconsistency Detection ===")
    
    # Find rows where fouling increases but U/Q also increase
    # Sort by R_f_hot and check if U/Q increase
    df_sorted_hot = df.sort_values('R_f_hot').reset_index(drop=True)
    inconsistencies_hot = []
    
    for i in range(1, len(df_sorted_hot)):
        prev = df_sorted_hot.iloc[i-1]
        curr = df_sorted_hot.iloc[i]
        
        if curr['R_f_hot'] > prev['R_f_hot']:
            if curr['U_overall'] > prev['U_overall']:
                inconsistencies_hot.append({
                    'type': 'R_f_hot_increase_U_increase',
                    'index': df_sorted_hot.index[i],
                    'R_f_hot': (prev['R_f_hot'], curr['R_f_hot']),
                    'U_overall': (prev['U_overall'], curr['U_overall']),
                })
            if curr['Q_actual'] > prev['Q_actual']:
                inconsistencies_hot.append({
                    'type': 'R_f_hot_increase_Q_increase',
                    'index': df_sorted_hot.index[i],
                    'R_f_hot': (prev['R_f_hot'], curr['R_f_hot']),
                    'Q_actual': (prev['Q_actual'], curr['Q_actual']),
                })
    
    # Same for R_f_cold
    df_sorted_cold = df.sort_values('R_f_cold').reset_index(drop=True)
    inconsistencies_cold = []
    
    for i in range(1, len(df_sorted_cold)):
        prev = df_sorted_cold.iloc[i-1]
        curr = df_sorted_cold.iloc[i]
        
        if curr['R_f_cold'] > prev['R_f_cold']:
            if curr['U_overall'] > prev['U_overall']:
                inconsistencies_cold.append({
                    'type': 'R_f_cold_increase_U_increase',
                    'index': df_sorted_cold.index[i],
                    'R_f_cold': (prev['R_f_cold'], curr['R_f_cold']),
                    'U_overall': (prev['U_overall'], curr['U_overall']),
                })
            if curr['Q_actual'] > prev['Q_actual']:
                inconsistencies_cold.append({
                    'type': 'R_f_cold_increase_Q_increase',
                    'index': df_sorted_cold.index[i],
                    'R_f_cold': (prev['R_f_cold'], curr['R_f_cold']),
                    'Q_actual': (prev['Q_actual'], curr['Q_actual']),
                })
    
    results['inconsistencies'] = inconsistencies_hot + inconsistencies_cold
    
    print(f"Found {len(inconsistencies_hot)} inconsistencies with R_f_hot")
    print(f"Found {len(inconsistencies_cold)} inconsistencies with R_f_cold")
    print(f"Total inconsistencies: {len(results['inconsistencies'])}")
    
    # Summary
    print("\n=== Summary ===")
    all_negative = all([
        corr_hot_U < 0,
        corr_hot_Q < 0,
        corr_cold_U < 0,
        corr_cold_Q < 0,
        corr_total_U < 0,
        corr_total_Q < 0,
    ])
    
    if all_negative and len(results['inconsistencies']) == 0:
        print("PASS: Dataset shows physically consistent behavior")
        print("All correlations are negative as expected.")
    else:
        print("WARNING: Dataset may have physics inconsistencies")
        if not all_negative:
            print("Some correlations are positive (should be negative)")
        if len(results['inconsistencies']) > 0:
            print(f"Found {len(results['inconsistencies'])} rows with inconsistent behavior")
    
    return results


def main():
    """Validate all dataset files in the data directory."""
    data_dir = 'data'
    
    if not os.path.exists(data_dir):
        print(f"Data directory not found: {data_dir}")
        return
    
    csv_files = [f for f in os.listdir(data_dir) if f.endswith('.csv')]
    
    if not csv_files:
        print(f"No CSV files found in {data_dir}")
        return
    
    print(f"Found {len(csv_files)} dataset file(s)\n")
    
    for csv_file in csv_files:
        csv_path = os.path.join(data_dir, csv_file)
        print("="*60)
        results = validate_dataset_physics(csv_path)
        print()


if __name__ == "__main__":
    main()

