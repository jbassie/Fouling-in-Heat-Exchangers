import random
import numpy as np
from config import TOLERANCE, MAX_ITER
from phe_model import PlateHeatExchanger
from fluids import FluidProperties


class Solver:
    def __init__(self):
        self.phe = PlateHeatExchanger()
        self.fluid_hot = FluidProperties('flue_gas')
        self.fluid_cold = FluidProperties('water')

    
    def run_simulation(self, inputs):
        """
        Executes the simulation for a given set of inputs using an iterative approach
        Inputs: Dictionary with flow rates, temperatures, and fouling factors
        """

        #Unpack Inputs(Step 1)
        m_hot = inputs['m_dot_hot']
        m_cold = inputs['m_dot_cold']
        T_hot_in = inputs['T_hot_in']
        T_cold_in = inputs['T_cold_in']
        Rf_hot = inputs['R_f_hot']
        Rf_cold = inputs['R_f_cold']


        #Step 2 : Initialize Guess (T_out_guess <= T_in)
        T_hot_out_guess = T_hot_in - 10 #initial guess for hot outlet temperature
        T_hot_out_curr = T_hot_in

        iteration = 0
        error = float('inf')

        #Step 3 : Iterative Solution
        while error > TOLERANCE and iteration < MAX_ITER:
            #3a Bulk Mean Temperatures(Hot)
            T_hot_mean = (T_hot_in +T_hot_out_guess ) / 2
            props_hot = self.fluid_hot.get_properties(T_hot_mean)

            #3b Hot Side Physics(Fouling & Hydrodynamics)
            res_h = self.phe._calculate_side_physics(m_hot, props_hot, Rf_hot)


           #3c Energy Balance for Cold Side Guess
           # Q_guess = m_h * Cp_h * (T_in - T_out_guess)
            # T_c_out = T_c_in + Q_guess / (m_c * Cp_c_approx)
            # We estimate Cp_c at inlet for the initial energy balance guess
            props_cold_in = self.fluid_cold.get_properties(T_cold_in)
            Q_approx = m_hot * props_hot['cp'] * (T_hot_in - T_hot_out_guess)
            T_cold_out_guess = T_cold_in + Q_approx / (m_cold * props_cold_in['cp'])

            T_c_mean = (T_cold_in + T_cold_out_guess) / 2
            props_cold = self.fluid_cold.get_properties(T_c_mean)

            #3d. Cold Side Physics(Fouling & Hydrodynamics)
            res_c = self.phe._calculate_side_physics(m_cold, props_cold, Rf_cold)

            #3e. Overall U, NTU, Epsillon Calculation
            U_overall = self.phe.calc_overall_U(res_h['h'], res_c['h'], Rf_hot, Rf_cold)
            epsilon, C_min, C_h, _ = self.phe.calc_epsilon_ntu(
                U_overall, m_hot, props_hot['cp'], m_cold, props_cold['cp']
                )
            
            #3f. Calculate Actual Outlet Temperatures
            Q_actual = epsilon * C_min * (T_hot_in - T_cold_in)
            T_hot_out_calc = T_hot_in - (Q_actual / C_h)

            #3g. Calculate Guess Error(Relaxation)
            T_h_out_guess = (T_hot_out_calc + T_hot_out_guess) / 2

            #check convergence
            error = abs(T_h_out_guess - T_hot_out_guess)
            iteration += 1

            T_h_out_final = T_h_out_calc
            #Calculate Water Outlet based on Converged Q

            T_c_out_final = T_cold_in + (Q_actual / (m_cold * props_cold['cp']))
        
        #Step 4 : Record Data
        return {
            'inputs': inputs,
            'results': {
                'T_hot_out': T_h_out_final,
                'T_cold_out': T_c_out_final,
                'Q_actual': Q_actual,
                'U_overall': U_overall,
                'h_hot': res_h['h'],
                'h_cold': res_c['h'],
                'error': error,
                 'iteration': iteration,
            }
        }