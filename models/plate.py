"""
Plate Heat Exchanger (PHE) Model
This file contains the model for the Plate Heat Exchanger (PHE).
This model is based on the literature and the geometry of the PHE.
"""
import numpy as np
from config import K_DEPOSIT, K_PLATE
from models.base_model import BaseHeatExchanger

class PlateHeatExchanger(BaseHeatExchanger):
    def __init__(self):
        #Geometry parameters (Fixed parameters for simulation)
        self.lengthOfPlate = 1.0               #Length of plate (m)
        self.widthofPlate = 0.5               #Width of plate (m)
        self.spacingBetweenPlates = 0.003     #Spacing/Gap between plates (m) - Very narrow!
        self.chevronAngle = 60                #Chevron angle (degrees)
        self.numberOfPlates = 50              #Total number of channels
        self.plateThickness = 0.0006          #Plate thickness (m)
        self.passes_hot = 1                   #Number of passes on hot side
        self.passes_cold = 1                  #Number of passes on cold side

        # --- Channel asymmetry (VERY IMPORTANT) ---
        # Cold side typically has more channels in heating applications
        # This is because the cold side is the side that is being heated
        self.channels_hot = int(0.45 * self.numberOfPlates)
        self.channels_cold = self.numberOfPlates - self.channels_hot

    def _calc_hydraulic_diameter(self, fouling_thickness):
        """
        Calculate the hydraulic diameter based on the reduced gap
        Adapts the logic of diameter reduction in literature (Eqns 2 -3)
        """
        b_eff = self.spacingBetweenPlates - (2 * fouling_thickness)

        if b_eff <= 0:
            b_eff = 1e-6  # Minimum effective gap to avoid division by zero
        return 2 * b_eff, b_eff

    def _calc_nusselt(self, Re, Pr):
        """
        Generic Chisholm & Wanniarachchi correlation for Chevron plates in PHE
        Replaces Zukauskas correlation(Eqn 13) and sieder-tate correlation(Eqn 7)
        """

        return 0.724 * (Re ** 0.583) * (Pr ** 0.333)
    
    def _calculate_side_physics(self, m_dot, props, Rf, side):
        """
        Computes velocity, Reynolds, and Heat Transfer Coefficient (h) 
        Corresponds to the 3b(Hot) and 3d(Cold) in Table 2
        
        Note: Uses unfouled geometry for flow area calculation to ensure fouling
        resistance always reduces heat transfer. Fouling is applied as thermal
        resistance in calc_overall_U(), not through geometry reduction.
        """
        #1. Calculate fouling layer thickness from Resistance (R = t/k)
        #Rearranged from Eq 4/11 logic (for reference/reporting only)
        t_fouling = Rf * K_DEPOSIT

        #2. Use unfouled geometry for flow area calculation
        # This prevents velocity increase from dominating fouling resistance
        # Fouling resistance is applied separately in calc_overall_U()
        b_eff_unfouled = self.spacingBetweenPlates
        Dh_eff = 2 * b_eff_unfouled

        # --- Channel count asymmetry ---
        if side == "hot":
            channels = self.channels_hot
            passes = self.passes_hot
            chevron_exp = 0.60   # slightly weaker mixing
        else:
            channels = self.channels_cold
            passes = self.passes_cold
            chevron_exp = 0.70   # stronger mixing on cold side

            # --- Pass-number effect ---
            # Fewer parallel channels per pass → higher velocity
        channels_per_pass = channels / passes

            # --- Flow area ---
        flow_area = self.widthofPlate * b_eff_unfouled * channels_per_pass

        #3 Hydrodynamics
        # In a PHE, channels alternate between hot and cold sides
       
        velocity  = m_dot / (props['rho'] * flow_area)
        Re = (props['rho'] * velocity * Dh_eff) / props['mu']

        chevron_factor = np.sin(np.radians(self.chevronAngle)) ** chevron_exp

        #4 Heat Transfer Coefficient
        Nu = chevron_factor * self._calc_nusselt(Re, props['pr'])
        h = (Nu * props['k']) / Dh_eff

        return {
            'h': h,
            'Velocity': velocity,
            'Re': Re,
            't_fouling': t_fouling,
        }
    
    def calc_overall_U(self, h_hot, h_cold, Rf_hot, Rf_cold):
        """
        Calculates overall heat transfer coefficient (U)
        Adapted Eq 14 for flat plates
        """

        R_wall = self.plateThickness / K_PLATE  #Thermal resistance of plate material
        R_total = (1/h_hot) + (1/h_cold) + R_wall +  Rf_hot + Rf_cold
        return 1/ R_total
    

    def calc_ntu(self, U, m_dot_hot, cp_hot, m_dot_cold, cp_cold):
        """
        Calculates NTU (Number of Transfer Units)
        Corresponds to the 4a in Table 2
        """
        C_h = m_dot_hot * cp_hot
        C_c = m_dot_cold * cp_cold
        C_min = min(C_h, C_c)

        Area_total = self.lengthOfPlate * self.widthofPlate * (self.numberOfPlates * 2)
        NTU = U * Area_total / C_min

        return NTU
    
    def calc_epsilon_ntu(self, U, m_dot_hot, cp_hot, m_dot_cold, cp_cold):
        """
        Calculates effectiveness (epsilon) and NTU (Number of Transfer Units)
        Corresponds to the 4a and 4b in Table 2
        
        Note: This method is kept for backward compatibility.
        Use calc_effectiveness() for the standard interface.
        """
        C_h = m_dot_hot * cp_hot
        C_c = m_dot_cold * cp_cold

        C_min = min(C_h, C_c)
        C_max = max(C_h, C_c)
        C_r = C_min / C_max

        Area_total = self.lengthOfPlate * self.widthofPlate * (self.numberOfPlates * 2)
        NTU = U * Area_total / C_min

        #Counter-Flow Effectiveness Equation
        if C_r < 1:
            numer = 1 - np.exp(-NTU * (1 - C_r))
            denom = 1 - C_r * np.exp(-NTU * (1 - C_r))
            epsilon = numer / denom
        else: # C_r = 1 case
            epsilon = NTU / (1 + NTU)

        return epsilon, C_min, C_h, C_c


    def calc_effectiveness(self, U, NTU, m_dot_hot, cp_hot, m_dot_cold, cp_cold):
        """
        Calculate heat exchanger effectiveness.
        Uses standard counter-flow effectiveness for PHE.
        
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
        C_h = m_dot_hot * cp_hot
        C_c = m_dot_cold * cp_cold
        C_min = min(C_h, C_c)
        C_max = max(C_h, C_c)
        C_r = C_min / C_max

        #Counter-Flow Effectiveness Equation
        if C_r < 1:
            numer = 1 - np.exp(-NTU * (1 - C_r))
            denom = 1 - C_r * np.exp(-NTU * (1 - C_r))
            epsilon = numer / denom
        else: # C_r = 1 case
            epsilon = NTU / (1 + NTU)

        return epsilon, C_min, C_h, C_c
    
    def calc_phe_effectiveness(
        self,
        U,
        m_dot_hot, cp_hot,
        m_dot_cold, cp_cold,
        Rf_hot, Rf_cold,
        chevron_angle,
        maldistribution_factor=0.85
    ):
        """
            PHE-specific effectiveness–NTU formulation accounting for:
            - Channel imbalance
            - Fouling contact resistance
            - Chevron angle degradation
            - Maldistribution
            - Plate conduction
            
        Note: This is a specialized method. For standard interface, use calc_effectiveness().
        """

        # --- Heat capacity rates ---
        C_h = m_dot_hot * cp_hot
        C_c = m_dot_cold * cp_cold
        C_min = min(C_h, C_c)
        C_max = max(C_h, C_c)
        C_r = C_min / C_max

        # --- Effective heat transfer area ---
        A_total = self.lengthOfPlate * self.widthofPlate * (self.numberOfPlates * 2)

        # --- Plate conduction resistance (non-negligible!) ---
        R_plate = self.plateThickness / K_PLATE

        # --- Fouling-induced contact degradation ---
        fouling_penalty = 1.0 / (1.0 + 5.0 * (Rf_hot + Rf_cold))

        # --- Chevron angle effectiveness factor ---
        # Optimal mixing around 60°
        chevron_factor = np.sin(np.radians(chevron_angle)) ** 0.7

        # --- Effective NTU ---
        NTU_eff = (
            U * A_total / C_min
            * chevron_factor
            * maldistribution_factor
            * fouling_penalty
            / (1.0 + R_plate * U)
        )

        # --- PHE mixed-flow effectiveness correlation ---
        # Literature-supported approximation for multi-pass PHEs
        if C_r < 1:
            epsilon = (
                1
                - np.exp(-NTU_eff * (1 - C_r))
            ) / (
                1
                - C_r * np.exp(-NTU_eff * (1 - C_r))
            )
        else:
            epsilon = NTU_eff / (1 + NTU_eff)

        return epsilon, NTU_eff, C_min, C_h, C_c