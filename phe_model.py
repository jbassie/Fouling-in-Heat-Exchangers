"""
Plate Heat Exchanger (PHE) Model
This file contains the model for the Plate Heat Exchanger (PHE).
This model is based on the literature and the geometry of the PHE.
"""
import numpy as np
from config import K_DEPOSIT, K_PLATE

class PlateHeatExchanger:
    def __init__(self):
        #Geometry parameters (Fixed parameters for simulation)
       self.lengthOfPlate = 1.0               #Length of plate (m)
       self.widthofPlate = 0.5               #Width of plate (m)
       self.spacingBetweenPlates = 0.003     #Spacing/Gap between plates (m) - Very narrow!
       self.chevronAngle = 60                  #Chevron angle (degrees)
       self.numberOfPlates= 50              #Total number of channels
       self.plateThickness = 0.0006                    #Plate thickness (m)

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
    
    def _calculate_side_physics(self, m_dot, props, Rf):
        """
        Computes velocity, Reynolds, and Heat Transfer Coefficient (h) 
        Corresponds to the 3b(Hot) and 3d(Cold) in Table 2
        """
        #1. Calculate fouling layer thickness from Resistance (R = t/k)
        #Rearranged from Eq 4/11 logic
        t_fouling = Rf * K_DEPOSIT

        #2. Geometry Adjustments(Fouling reduces flow area)
        Dh_eff, b_eff = self._calc_hydraulic_diameter(t_fouling)

        #3 Hydrodynamics
        flow_area = self.widthofPlate * b_eff * self.numberOfPlates
        velocity  = m_dot / (props['rho'] * flow_area)
        Re = (props['rho'] * velocity * Dh_eff) / props['mu']

        #4 Heat Transfer Coefficient
        Nu = self._calc_nusselt(Re, props['pr'])
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
    

    def calc_epsilon_ntu(self, U, m_dot_hot, cp_hot, m_dot_cold, cp_cold):
        """
        Calculates effectiveness (epsilon) and NTU (Number of Transfer Units)
        Corresponds to the 4a and 4b in Table 2
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
