"""
Base Heat Exchanger Class

Provides common functionality shared between different heat exchanger types.
Implements DRY principle by centralizing shared methods.
"""
from abc import ABC, abstractmethod
from config import K_DEPOSIT, K_PLATE


class BaseHeatExchanger(ABC):
    """
    Abstract base class for heat exchanger models.
    
    Provides common methods and defines interface for exchanger-specific calculations.
    """
    
    def _calc_fouling_thickness(self, Rf: float) -> float:
        """
        Calculate fouling layer thickness from resistance.
        Equation: R = t/k, so t = R * k
        
        Args:
            Rf: Fouling resistance (m²·K/W)
            
        Returns:
            Fouling thickness (m)
        """
        return Rf * K_DEPOSIT
    
    @abstractmethod
    def calc_overall_U(self, h_hot: float, h_cold: float, Rf_hot: float, Rf_cold: float) -> float:
        """
        Calculate overall heat transfer coefficient.
        
        Args:
            h_hot: Hot side heat transfer coefficient (W/(m²·K))
            h_cold: Cold side heat transfer coefficient (W/(m²·K))
            Rf_hot: Hot side fouling resistance (m²·K/W)
            Rf_cold: Cold side fouling resistance (m²·K/W)
            
        Returns:
            Overall heat transfer coefficient (W/(m²·K))
        """
        pass
    
    @abstractmethod
    def calc_ntu(self, U: float, m_dot_hot: float, cp_hot: float,
                 m_dot_cold: float, cp_cold: float) -> float:
        """
        Calculate Number of Transfer Units (NTU).
        
        Args:
            U: Overall heat transfer coefficient (W/(m²·K))
            m_dot_hot: Hot side mass flow rate (kg/s)
            cp_hot: Hot side specific heat (J/(kg·K))
            m_dot_cold: Cold side mass flow rate (kg/s)
            cp_cold: Cold side specific heat (J/(kg·K))
            
        Returns:
            NTU (dimensionless)
        """
        pass
    
    @abstractmethod
    def calc_effectiveness(self, U: float, NTU: float, m_dot_hot: float, cp_hot: float,
                          m_dot_cold: float, cp_cold: float) -> tuple:
        """
        Calculate heat exchanger effectiveness.
        
        Args:
            U: Overall heat transfer coefficient (W/(m²·K))
            NTU: Number of Transfer Units
            m_dot_hot: Hot side mass flow rate (kg/s)
            cp_hot: Hot side specific heat (J/(kg·K))
            m_dot_cold: Cold side mass flow rate (kg/s)
            cp_cold: Cold side specific heat (J/(kg·K))
            
        Returns:
            Tuple of (effectiveness, C_min, C_h, C_c) where:
            - effectiveness: dimensionless (0-1)
            - C_min: minimum heat capacity rate (W/K)
            - C_h: hot side heat capacity rate (W/K)
            - C_c: cold side heat capacity rate (W/K)
        """
        pass

