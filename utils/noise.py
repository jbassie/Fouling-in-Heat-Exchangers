"""
Noise Generation Utilities

Functions for adding realistic noise to simulation data.
"""
from typing import Optional, List
import numpy as np
import pandas as pd
from config import NOISE_LEVEL, PERCENTAGE_NOISE, NOISE_COLUMNS


def add_noise_to_dataframe(df: pd.DataFrame, noise_level: Optional[float] = None,
                          percentage_noise: Optional[float] = None,
                          noise_columns: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Add Gaussian noise to a random percentage of rows in specified columns.
    
    Args:
        df: Input dataframe
        noise_level: Standard deviation as fraction of mean (default from config NOISE_LEVEL)
        percentage_noise: Percentage of rows to add noise to (default from config PERCENTAGE_NOISE)
        noise_columns: List of column names to add noise to (default from config)
        
    Returns:
        Dataframe with noise added to specified columns for selected rows
    """
    if noise_level is None:
        noise_level = NOISE_LEVEL
    if percentage_noise is None:
        percentage_noise = PERCENTAGE_NOISE
    if noise_columns is None:
        noise_columns = NOISE_COLUMNS
    
    df_noisy = df.copy()
    
    # Calculate number of rows to add noise to
    total_rows = len(df_noisy)
    num_noisy_rows = int(total_rows * percentage_noise)
    
    # Randomly select rows to add noise to
    np.random.seed()  # Use random seed for each call
    noisy_row_indices = np.random.choice(total_rows, size=num_noisy_rows, replace=False)
    
    for col in noise_columns:
        if col in df_noisy.columns:
            # Calculate standard deviation as fraction of mean value
            mean_val = df_noisy[col].abs().mean()
            if mean_val > 0:
                std_dev = mean_val * noise_level
                # Add Gaussian noise only to selected rows
                noise = np.zeros(len(df_noisy))
                noise[noisy_row_indices] = np.random.normal(0, std_dev, size=num_noisy_rows)
                df_noisy[col] = df_noisy[col] + noise
    
    return df_noisy


class NoiseInjector:
    """
    Adds realistic, traceable noise to simulation data.
    """

    def __init__(
        self,
        noise_columns: List[str],
        base_noise_frac: float = 0.002,
        spike_noise_frac: float = 0.01,
        spike_row_fraction: float = 0.1,
        random_state: Optional[int] = None,
    ):
        """
        Args:
            noise_columns: columns to inject noise into
            base_noise_frac: small noise applied to all rows (fraction of mean)
            spike_noise_frac: stronger noise applied to subset of rows
            spike_row_fraction: fraction of rows receiving spike noise
            random_state: reproducibility
        """
        self.noise_columns = noise_columns
        self.base_noise_frac = base_noise_frac
        self.spike_noise_frac = spike_noise_frac
        self.spike_row_fraction = spike_row_fraction
        self.rng = np.random.default_rng(random_state)

    def apply(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        n = len(df)

        # --- select spike rows ---
        num_spike = int(n * self.spike_row_fraction)
        spike_idx = self.rng.choice(n, size=num_spike, replace=False)

        df["noise_applied"] = 0
        df.loc[spike_idx, "noise_applied"] = 1

        for col in self.noise_columns:
            if col not in df.columns:
                continue

            mean_val = df[col].abs().mean()
            if mean_val == 0:
                continue

            # --- noise scales ---
            base_std = mean_val * self.base_noise_frac
            spike_std = mean_val * self.spike_noise_frac

            # --- generate noise ---
            base_noise = self.rng.normal(0, base_std, size=n)
            spike_noise = np.zeros(n)
            spike_vals = self.rng.normal(0, spike_std, size=num_spike)
            spike_noise[spike_idx] = spike_vals

            total_noise = base_noise + spike_noise

            # --- apply ---
            df[col] += total_noise

            # --- tracking ---
            df[f"noise_{col}"] = 0
            df.loc[spike_idx, f"noise_{col}"] = 1

            df[f"delta_{col}"] = total_noise

        return df
