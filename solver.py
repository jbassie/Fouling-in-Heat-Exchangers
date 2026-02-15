"""
Unified Solver for Heat Exchanger Simulations

Works with both Plate Heat Exchanger (PHE) and Cross-Flow Heat Exchanger types.
Implements the iterative solution algorithm common to both exchanger types.
"""
from config import TOLERANCE, MAX_ITER
from fluids import FluidProperties
from models.base_model import BaseHeatExchanger


class Solver:
    """
    Unified solver that works with any heat exchanger type.
    
    Implements the iterative solution algorithm:
    1. Random selection of inputs (done externally)
    2. Initialize guess for outlet temperatures
    3. Iterative solution until convergence
    4. Record data with fouling resistances
    """
    
    def __init__(self, exchanger: BaseHeatExchanger):
        """
        Initialize solver with a heat exchanger instance.
        
        Args:
            exchanger: Instance of BaseHeatExchanger (PHE or CrossFlow)
        """
        self.exchanger = exchanger
        self.fluid_hot = FluidProperties('flue_gas')
        self.fluid_cold = FluidProperties('water')
    
    def _calc_hot_side_physics(self, m_dot: float, props: dict, Rf: float) -> dict:
        """
        Calculate hot side (flue gas) heat transfer coefficient and related variables.
        
        Delegates to exchanger-specific methods.
        
        Args:
            m_dot: Mass flow rate (kg/s)
            props: Fluid properties dictionary
            Rf: Fouling resistance (m²·K/W)
            
        Returns:
            Dictionary with 'h' and other side-specific variables
        """
        # Check if exchanger has specific method for hot side
        if hasattr(self.exchanger, '_calc_h_flue'):
            # CrossFlow exchanger
            return self.exchanger._calc_h_flue(m_dot, props, Rf)
        elif hasattr(self.exchanger, '_calculate_side_physics'):
            # PHE exchanger
            return self.exchanger._calculate_side_physics(m_dot, props, Rf, 'hot')
        else:
            raise NotImplementedError("Exchanger must implement hot side physics calculation")
    
    def _calc_cold_side_physics(self, m_dot: float, props: dict, Rf: float) -> dict:
        """
        Calculate cold side (water) heat transfer coefficient and related variables.
        
        Delegates to exchanger-specific methods.
        
        Args:
            m_dot: Mass flow rate (kg/s)
            props: Fluid properties dictionary
            Rf: Fouling resistance (m²·K/W)
            
        Returns:
            Dictionary with 'h' and other side-specific variables
        """
        # Check if exchanger has specific method for cold side
        if hasattr(self.exchanger, '_calc_h_water'):
            # CrossFlow exchanger
            return self.exchanger._calc_h_water(m_dot, props, Rf)
        elif hasattr(self.exchanger, '_calculate_side_physics'):
            # PHE exchanger
            return self.exchanger._calculate_side_physics(m_dot, props, Rf, 'cold')
        else:
            raise NotImplementedError("Exchanger must implement cold side physics calculation")
    
    def run_simulation(self, inputs: dict) -> dict:
        """
        Run heat exchanger simulation following the iterative algorithm.
        
        Algorithm:
        (1) Randomly select inlet flow-rates, temperatures, and fouling factors
        (2) Guess T_fo,guess ≤ T_fi and initialize T_fo = T_fi
        (3) Repeat (a)-(g) until |T_fo,guess - T_fo| < 10^-6:
            (a) Evaluate flue gas properties at bulk mean temperature
            (b) Compute flue-side fouling variables (D_f,o, V_max,f, M_dot_flue,f, h_flue)
            (c) Evaluate T_wo,guess from energy balance and water properties at bulk mean
            (d) Determine water-side fouling variables (D_f,i, V_f, M_dot_w,f, h_w)
            (e) Calculate U = f(h_flue, h_w, R1, R2), NTU = f(U, (M*Cp)_min), ε = f(U, NTU, Cr)
            (f) Obtain T_fo and T_wo from computations
            (g) Update T_fo,guess = (T_fo,guess + T_fo)/2 and T_wo,guess = (T_wo,guess + T_wo)/2
        (4) Record generated data with fouling resistances (R1, R2)
        
        Args:
            inputs: Dictionary with keys:
                - 'm_dot_hot': Hot side mass flow rate (kg/s)
                - 'm_dot_cold': Cold side mass flow rate (kg/s)
                - 'T_hot_in': Hot side inlet temperature (°C)
                - 'T_cold_in': Cold side inlet temperature (°C)
                - 'R_f_hot': Hot side fouling resistance (m²·K/W)
                - 'R_f_cold': Cold side fouling resistance (m²·K/W)
        
        Returns:
            Dictionary with inputs and results
        """
        # Step 1: Unpack inputs
        m_dot_hot = inputs['m_dot_hot']
        m_dot_cold = inputs['m_dot_cold']
        T_hot_in = inputs['T_hot_in']
        T_cold_in = inputs['T_cold_in']
        Rf_hot = inputs['R_f_hot']
        Rf_cold = inputs['R_f_cold']
        
        # Step 2: Initialize guess
        # Algorithm: T_fo,guess ≤ T_fi and initialize T_fo = T_fi
        T_hot_out_guess = T_hot_in - 10.0  # Initial guess (T_out <= T_in)
        T_hot_out = T_hot_in  # Initialize T_fo = T_fi (as per algorithm)
        T_cold_out_guess = T_cold_in  # Initialize cold side guess
        
        iteration = 0
        error = float('inf')
        
        # Step 3: Iterative solution
        while error > TOLERANCE and iteration < MAX_ITER:
            # (a) Evaluate hot side properties at bulk mean temperature
            T_hot_mean = (T_hot_in + T_hot_out_guess) / 2
            props_hot = self.fluid_hot.get_properties(T_hot_mean)
            
            # (b) Compute hot-side fouling variables and h_hot
            res_hot = self._calc_hot_side_physics(m_dot_hot, props_hot, Rf_hot)
            h_hot = res_hot['h']
            
            # (c) Evaluate cold side outlet guess from energy balance
            # Energy balance: Q = m_hot * cp_hot * (T_hot_in - T_hot_out_guess)
            Q_guess = m_dot_hot * props_hot['cp'] * (T_hot_in - T_hot_out_guess)
            
            # Estimate cold outlet using inlet properties
            props_cold_in = self.fluid_cold.get_properties(T_cold_in)
            T_cold_out_guess = T_cold_in + Q_guess / (m_dot_cold * props_cold_in['cp'])
            
            # Cold side properties at bulk mean temperature
            T_cold_mean = (T_cold_in + T_cold_out_guess) / 2
            props_cold = self.fluid_cold.get_properties(T_cold_mean)
            
            # (d) Determine cold-side fouling variables and h_cold
            res_cold = self._calc_cold_side_physics(m_dot_cold, props_cold, Rf_cold)
            h_cold = res_cold['h']
            
            # (e) Calculate U, NTU, epsilon
            U = self.exchanger.calc_overall_U(h_hot, h_cold, Rf_hot, Rf_cold)
            NTU = self.exchanger.calc_ntu(U, m_dot_hot, props_hot['cp'],
                                         m_dot_cold, props_cold['cp'])
            
            # Use PHE-specific effectiveness if available, otherwise standard
            if hasattr(self.exchanger, 'calc_phe_effectiveness'):
                chevron_angle = getattr(self.exchanger, 'chevronAngle', 60)
                epsilon, NTU_eff, C_min, C_h, C_c = self.exchanger.calc_phe_effectiveness(
                    U, m_dot_hot, props_hot['cp'],
                    m_dot_cold, props_cold['cp'],
                    Rf_hot, Rf_cold,
                    chevron_angle
                )
                NTU = NTU_eff  # Use effective NTU for PHE
            else:
                epsilon, C_min, C_h, C_c = self.exchanger.calc_effectiveness(
                    U, NTU, m_dot_hot, props_hot['cp'],
                    m_dot_cold, props_cold['cp']
                )
            
            # (f) Obtain actual outlet temperatures
            Q_actual = epsilon * C_min * (T_hot_in - T_cold_in)
            T_hot_out_calc = T_hot_in - (Q_actual / C_h)
            T_cold_out_calc = T_cold_in + (Q_actual / C_c)
            
            # (g) Update guesses using relaxation
            # Algorithm: T_fo,guess = (T_fo,guess + T_fo)/2 and T_wo,guess = (T_wo,guess + T_wo)/2
            T_hot_out_guess = (T_hot_out_guess + T_hot_out_calc) / 2
            T_cold_out_guess = (T_cold_out_guess + T_cold_out_calc) / 2
            T_hot_out = T_hot_out_calc  # Update T_fo
            
            # Check convergence: |T_fo,guess - T_fo| < 10^-6
            error = abs(T_hot_out_guess - T_hot_out)
            iteration += 1
        
        # Step 4: Record data
        # Use final calculated values from convergence
        # Recalculate cold side with final temperatures for accuracy (using final Q_actual)
        T_cold_mean_final = (T_cold_in + T_cold_out_guess) / 2
        props_cold_final = self.fluid_cold.get_properties(T_cold_mean_final)
        T_cold_out_final = T_cold_in + (Q_actual / (m_dot_cold * props_cold_final['cp']))
        
        # Build results dictionary
        results = {
            'T_hot_out': T_hot_out,  # Use T_fo (final converged value)
            'T_cold_out': T_cold_out_final,
            'Q_actual': Q_actual,
            'U_overall': U,
            'NTU': NTU,
            'epsilon': epsilon,
            'h_hot': h_hot,
            'h_cold': h_cold,
            'error': error,
            'iteration': iteration,
        }
        
        # Add exchanger-specific results
        if 'D_f_o' in res_hot:
            results['D_f_o'] = res_hot['D_f_o']
        if 'D_f_i' in res_cold:
            results['D_f_i'] = res_cold['D_f_i']
        if 'V_max_f' in res_hot:
            results['V_max_f'] = res_hot['V_max_f']
        if 'V_f' in res_cold:
            results['V_f'] = res_cold['V_f']
        if 'M_dot_f' in res_hot:
            results['M_dot_flue_f'] = res_hot['M_dot_f']
        if 'M_dot_f' in res_cold:
            results['M_dot_water_f'] = res_cold['M_dot_f']
        if 'Velocity' in res_hot:
            results['Velocity_hot'] = res_hot['Velocity']
        if 'Velocity' in res_cold:
            results['Velocity_cold'] = res_cold['Velocity']
        
        return {
            'inputs': inputs,
            'results': results
        }
