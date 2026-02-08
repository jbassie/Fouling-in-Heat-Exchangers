"""
Noise Generation Utilities

Functions for adding realistic noise to simulation data.
"""
import numpy as np
import pandas as pd
from config import NOISE_LEVEL, NOISE_COLUMNS


def add_noise_to_dataframe(df: pd.DataFrame, noise_level: float = None, 
                          noise_columns: list = None) -> pd.DataFrame:
    """
    Add Gaussian noise to specified columns in a dataframe.
    
    Args:
        df: Input dataframe
        noise_level: Standard deviation as fraction of mean (default from config)
        noise_columns: List of column names to add noise to (default from config)
        
    Returns:
        Dataframe with noise added to specified columns
    """
    if noise_level is None:
        noise_level = NOISE_LEVEL
    if noise_columns is None:
        noise_columns = NOISE_COLUMNS
    
    df_noisy = df.copy()
    
    for col in noise_columns:
        if col in df_noisy.columns:
            # Calculate standard deviation as fraction of mean value
            mean_val = df_noisy[col].abs().mean()
            if mean_val > 0:
                std_dev = mean_val * noise_level
                # Add Gaussian noise
                noise = np.random.normal(0, std_dev, size=len(df_noisy))
                df_noisy[col] = df_noisy[col] + noise
    
    return df_noisy

