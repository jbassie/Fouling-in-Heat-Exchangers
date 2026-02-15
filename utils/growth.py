"""
Fouling Growth Models

Supports:
- Linear growth
- Exponential growth
- Asymptotic (saturation-limited) growth

Designed for hourly time-stepping in heat exchanger simulations.
"""

import random
import numpy as np


def fouling_growth_model(
    current_Rf: float,
    growth_type: str = "linear",
    growth_rate_range=(1e-6, 5e-6),
    dt: float = 1.0,
    Rf_max: float = 0.001,
    tau: float = 1000.0,
):
    """
    Update fouling resistance based on selected growth model.

    Parameters
    ----------
    current_Rf : float
        Current fouling resistance (m²·K/W)

    growth_type : str
        'linear', 'exponential', or 'asymptotic'

    growth_rate_range : tuple
        Min/max fouling growth rate per hour

    dt : float
        Time step (hours)

    Rf_max : float
        Maximum fouling resistance (used for asymptotic model)

    tau : float
        Fouling time constant (hours) for asymptotic model

    Returns
    -------
    new_Rf : float
        Updated fouling resistance
    """

    growth_rate = random.uniform(*growth_rate_range)

    if growth_type.lower() == "linear":
        # Rf = Rf + k·dt
        new_Rf = current_Rf + growth_rate * dt

    elif growth_type.lower() == "exponential":
        # Rf = Rf · exp(k·dt)
        new_Rf = current_Rf * np.exp(growth_rate * dt)

    elif growth_type.lower() == "asymptotic":
        # Rf = Rf + (Rf_max - Rf) · (1 - exp(-dt / tau))
        #(Kern–Seaton–style) 
        new_Rf = current_Rf + (Rf_max - current_Rf) * (1 - np.exp(-dt / tau))

    else:
        raise ValueError(
            "Invalid growth_type. Choose from: 'linear', 'exponential', 'asymptotic'"
        )

    # Safety clamp (never exceed physical maximum)
    return min(new_Rf, Rf_max)
