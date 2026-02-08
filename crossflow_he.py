"""
Cross-Flow Heat Exchanger Model

This class implements a cross-flow heat exchanger model following the algorithm:
1. Random selection of inputs (done externally)
2. Iterative solution for outlet temperatures
3. Calculation of fouling variables and heat transfer coefficients
4. Overall U, NTU, and effectiveness calculation
"""
import numpy as np
from config import K_DEPOSIT, K_PLATE, TOLERANCE, MAX_ITER
from fluids import FluidProperties


class CrossFlowHeatExchanger:
    """
    Cross-flow heat exchanger model with fouling effects.
    
    Implements the algorithm:
    - Equations (4), (5), (6), (7) for water side (cold side)
    - Equations (10), (11), (12), (13) for flue gas side (hot side)
    - Cross-flow effectiveness-NTU method
    """
    
    def __init__(self):
        """Initialize cross-flow heat exchanger geometry parameters."""
        # Geometry parameters for cross-flow heat exchanger
        # Typical shell-and-tube or finned-tube configuration
        self.tube_outer_diameter = 0.025  # Outer diameter of tubes (m)
        self.tube_inner_diameter = 0.020  # Inner diameter of tubes (m)
        self.tube_length = 2.0  # Length of tubes (m)
        self.number_of_tubes = 100  # Number of tubes
        self.shell_diameter = 0.5  # Shell diameter (m)
        
        # Initialize fluid property handlers
        self.fluid_hot = FluidProperties('flue_gas')
        self.fluid_cold = FluidProperties('water')
    
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
    
    def _calc_hydraulic_diameter_water(self, D_f_i: float) -> float:
        """
        Calculate effective hydraulic diameter for water side (inside tubes).
        Equation (4): D_f,i - fouling reduces inner diameter
        
        Args:
            D_f_i: Fouling thickness on inner surface (m)
            
        Returns:
            Effective hydraulic diameter (m)
        """
        D_i_eff = self.tube_inner_diameter - (2 * D_f_i)
        if D_i_eff <= 0:
            D_i_eff = 1e-6  # Minimum to avoid division by zero
        return D_i_eff
    
    def _calc_velocity_water(self, m_dot: float, rho: float, D_f_i: float) -> float:
        """
        Calculate water velocity in tubes.
        Equation (5): V_f - velocity with fouling
        
        Args:
            m_dot: Mass flow rate (kg/s)
            rho: Density (kg/m³)
            D_f_i: Fouling thickness on inner surface (m)
            
        Returns:
            Velocity (m/s)
        """
        D_i_eff = self._calc_hydraulic_diameter_water(D_f_i)
        A_cross = np.pi * (D_i_eff / 2) ** 2
        A_total = A_cross * self.number_of_tubes
        return m_dot / (rho * A_total)
    
    def _calc_mass_flux_water(self, m_dot: float, D_f_i: float) -> float:
        """
        Calculate mass flux for water side.
        Equation (6): M_dot_w,f - mass flux with fouling
        
        Args:
            m_dot: Mass flow rate (kg/s)
            D_f_i: Fouling thickness on inner surface (m)
            
        Returns:
            Mass flux (kg/(m²·s))
        """
        D_i_eff = self._calc_hydraulic_diameter_water(D_f_i)
        A_cross = np.pi * (D_i_eff / 2) ** 2
        A_total = A_cross * self.number_of_tubes
        return m_dot / A_total
    
    def _calc_h_water(self, m_dot: float, props: dict, Rf: float) -> dict:
        """
        Calculate water-side heat transfer coefficient.
        Equations (4), (5), (6), (7) - water side fouling variables
        
        Args:
            m_dot: Mass flow rate (kg/s)
            props: Fluid properties dictionary
            Rf: Fouling resistance (m²·K/W)
            
        Returns:
            Dictionary with h, D_f_i, V_f, M_dot_f, and other variables
        """
        # Equation (4): D_f,i - fouling thickness
        D_f_i = self._calc_fouling_thickness(Rf)
        
        # Equation (5): V_f - velocity
        V_f = self._calc_velocity_water(m_dot, props['rho'], D_f_i)
        
        # Equation (6): M_dot_w,f - mass flux
        M_dot_f = self._calc_mass_flux_water(m_dot, D_f_i)
        
        # Equation (7): h_w - heat transfer coefficient
        # Use Dittus-Boelter correlation for turbulent flow in tubes
        D_i_eff = self._calc_hydraulic_diameter_water(D_f_i)
        Re = (props['rho'] * V_f * D_i_eff) / props['mu']
        
        if Re > 2300:  # Turbulent flow
            Nu = 0.023 * (Re ** 0.8) * (props['pr'] ** 0.4)
        else:  # Laminar flow
            Nu = 3.66  # Constant Nusselt for fully developed laminar flow
        
        h = (Nu * props['k']) / D_i_eff
        
        return {
            'h': h,
            'D_f_i': D_f_i,
            'V_f': V_f,
            'M_dot_f': M_dot_f,
            'Re': Re,
            'Nu': Nu,
        }
    
    def _calc_hydraulic_diameter_flue(self, D_f_o: float) -> float:
        """
        Calculate effective hydraulic diameter for flue gas side (outside tubes).
        Equation (10): D_f,o - fouling on outer surface
        
        Args:
            D_f_o: Fouling thickness on outer surface (m)
            
        Returns:
            Effective hydraulic diameter (m)
        """
        # For shell side, use equivalent diameter based on tube arrangement
        # Simplified: use shell diameter minus tube blockage
        D_o_eff = self.tube_outer_diameter + (2 * D_f_o)
        # Shell-side equivalent diameter (simplified)
        D_eq = self.shell_diameter - (self.number_of_tubes * D_o_eff / 2)
        if D_eq <= 0:
            D_eq = 1e-6
        return D_eq
    
    def _calc_velocity_max_flue(self, m_dot: float, rho: float, D_f_o: float) -> float:
        """
        Calculate maximum flue gas velocity in shell.
        Equation (11): V_max,f - maximum velocity with fouling
        
        Args:
            m_dot: Mass flow rate (kg/s)
            rho: Density (kg/m³)
            D_f_o: Fouling thickness on outer surface (m)
            
        Returns:
            Maximum velocity (m/s)
        """
        # Shell-side flow area (simplified)
        D_o_eff = self.tube_outer_diameter + (2 * D_f_o)
        A_shell = np.pi * (self.shell_diameter / 2) ** 2
        A_tubes = self.number_of_tubes * np.pi * (D_o_eff / 2) ** 2
        A_flow = A_shell - A_tubes
        
        if A_flow <= 0:
            A_flow = 1e-6
        
        return m_dot / (rho * A_flow)
    
    def _calc_mass_flux_flue(self, m_dot: float, D_f_o: float) -> float:
        """
        Calculate mass flux for flue gas side.
        Equation (12): M_dot_flue,f - mass flux with fouling
        
        Args:
            m_dot: Mass flow rate (kg/s)
            D_f_o: Fouling thickness on outer surface (m)
            
        Returns:
            Mass flux (kg/(m²·s))
        """
        D_o_eff = self.tube_outer_diameter + (2 * D_f_o)
        A_shell = np.pi * (self.shell_diameter / 2) ** 2
        A_tubes = self.number_of_tubes * np.pi * (D_o_eff / 2) ** 2
        A_flow = A_shell - A_tubes
        
        if A_flow <= 0:
            A_flow = 1e-6
        
        return m_dot / A_flow
    
    def _calc_h_flue(self, m_dot: float, props: dict, Rf: float) -> dict:
        """
        Calculate flue gas side heat transfer coefficient.
        Equations (10), (11), (12), (13) - flue gas side fouling variables
        
        Args:
            m_dot: Mass flow rate (kg/s)
            props: Fluid properties dictionary
            Rf: Fouling resistance (m²·K/W)
            
        Returns:
            Dictionary with h, D_f_o, V_max_f, M_dot_f, and other variables
        """
        # Equation (10): D_f,o - fouling thickness
        D_f_o = self._calc_fouling_thickness(Rf)
        
        # Equation (11): V_max,f - maximum velocity
        V_max_f = self._calc_velocity_max_flue(m_dot, props['rho'], D_f_o)
        
        # Equation (12): M_dot_flue,f - mass flux
        M_dot_f = self._calc_mass_flux_flue(m_dot, D_f_o)
        
        # Equation (13): h_flue - heat transfer coefficient
        # Use Zukauskas correlation for flow across tube banks
        D_o_eff = self.tube_outer_diameter + (2 * D_f_o)
        Re = (props['rho'] * V_max_f * D_o_eff) / props['mu']
        
        # Zukauskas correlation for staggered tube banks
        if Re < 1000:
            C = 0.9
            m = 0.4
        elif Re < 200000:
            C = 0.27
            m = 0.63
        else:
            C = 0.021
            m = 0.84
        
        Nu = C * (Re ** m) * (props['pr'] ** 0.36)
        h = (Nu * props['k']) / D_o_eff
        
        return {
            'h': h,
            'D_f_o': D_f_o,
            'V_max_f': V_max_f,
            'M_dot_f': M_dot_f,
            'Re': Re,
            'Nu': Nu,
        }
    
    def calc_overall_U(self, h_flue: float, h_water: float, R1: float, R2: float) -> float:
        """
        Calculate overall heat transfer coefficient.
        U = f(h_flue, h_water, R1, R2)
        
        Args:
            h_flue: Flue gas side heat transfer coefficient (W/(m²·K))
            h_water: Water side heat transfer coefficient (W/(m²·K))
            R1: Flue gas side fouling resistance (m²·K/W)
            R2: Water side fouling resistance (m²·K/W)
            
        Returns:
            Overall heat transfer coefficient (W/(m²·K))
        """
        # Wall resistance (tube wall)
        R_wall = (self.tube_outer_diameter - self.tube_inner_diameter) / (2 * K_PLATE)
        
        # Overall resistance (accounting for inner/outer area ratio)
        A_o = np.pi * self.tube_outer_diameter * self.tube_length * self.number_of_tubes
        A_i = np.pi * self.tube_inner_diameter * self.tube_length * self.number_of_tubes
        area_ratio = A_o / A_i
        
        R_total = (1 / h_flue) + R1 + R_wall + (area_ratio / h_water) + (area_ratio * R2)
        
        return 1 / R_total
    
    def calc_ntu(self, U: float, m_dot_hot: float, cp_hot: float, 
                  m_dot_cold: float, cp_cold: float) -> float:
        """
        Calculate Number of Transfer Units (NTU).
        NTU = f(U, (M*Cp)_min)
        
        Args:
            U: Overall heat transfer coefficient (W/(m²·K))
            m_dot_hot: Hot side mass flow rate (kg/s)
            cp_hot: Hot side specific heat (J/(kg·K))
            m_dot_cold: Cold side mass flow rate (kg/s)
            cp_cold: Cold side specific heat (J/(kg·K))
            
        Returns:
            NTU (dimensionless)
        """
        C_h = m_dot_hot * cp_hot
        C_c = m_dot_cold * cp_cold
        C_min = min(C_h, C_c)
        
        # Heat transfer area (outer surface area)
        A = np.pi * self.tube_outer_diameter * self.tube_length * self.number_of_tubes
        
        NTU = (U * A) / C_min
        
        return NTU
    
    def calc_effectiveness(self, U: float, NTU: float, m_dot_hot: float, cp_hot: float,
                          m_dot_cold: float, cp_cold: float) -> float:
        """
        Calculate heat exchanger effectiveness.
        epsilon = f(U, NTU, C_r)
        
        Uses cross-flow effectiveness correlation (both fluids unmixed).
        
        Args:
            U: Overall heat transfer coefficient (W/(m²·K))
            NTU: Number of Transfer Units
            m_dot_hot: Hot side mass flow rate (kg/s)
            cp_hot: Hot side specific heat (J/(kg·K))
            m_dot_cold: Cold side mass flow rate (kg/s)
            cp_cold: Cold side specific heat (J/(kg·K))
            
        Returns:
            Effectiveness (dimensionless, 0-1)
        """
        C_h = m_dot_hot * cp_hot
        C_c = m_dot_cold * cp_cold
        C_min = min(C_h, C_c)
        C_max = max(C_h, C_c)
        C_r = C_min / C_max
        
        # Cross-flow effectiveness (both fluids unmixed)
        # Approximation: epsilon = 1 - exp(-NTU^0.22 * (exp(-C_r * NTU^0.78) - 1) / C_r)
        if C_r < 1:
            epsilon = 1 - np.exp(
                (-NTU ** 0.22) * (np.exp(-C_r * NTU ** 0.78) - 1) / C_r
            )
        else:  # C_r = 1
            epsilon = 1 - np.exp(-NTU)
        
        return epsilon
    
    def run_simulation(self, inputs: dict) -> dict:
        """
        Run cross-flow heat exchanger simulation following the specified algorithm.
        
        Algorithm:
        (1) Random selection of inputs (done externally)
        (2) Guess T_fo,guess <= T_fi and initialize T_fo = T_fi
        (3) Repeat (a)-(g) till |T_fo,guess - T_fo| < 1e-6
        (a) Evaluate flue gas properties at bulk mean temperature
        (b) Compute flue-side fouling variables D_f,o, V_max,f, M_flue,f, h_flue (Eqns 10-13)
        (c) Evaluate T_wo,guess from energy balance and water properties at bulk mean
        (d) Determine water-side fouling variables D_f,i, V_f, M_w,f, h_w (Eqns 4-7)
        (e) Calculate U, NTU, epsilon
        (f) Obtain T_fo and T_wo from computations
        (g) Update T_fo,guess = (T_fo,guess + T_fo)/2 and T_wo,guess = (T_wo,guess + T_wo)/2
        (4) Record data with fouling resistances R1, R2
        
        Args:
            inputs: Dictionary with keys:
                - 'm_dot_hot': Flue gas mass flow rate (kg/s)
                - 'm_dot_cold': Water mass flow rate (kg/s)
                - 'T_hot_in': Flue gas inlet temperature (T_fi) (°C)
                - 'T_cold_in': Water inlet temperature (T_wi) (°C)
                - 'R_f_hot': Flue gas side fouling resistance (R1) (m²·K/W)
                - 'R_f_cold': Water side fouling resistance (R2) (m²·K/W)
        
        Returns:
            Dictionary with inputs and results
        """
        # Step 1: Unpack inputs
        m_dot_hot = inputs['m_dot_hot']
        m_dot_cold = inputs['m_dot_cold']
        T_fi = inputs['T_hot_in']  # Flue gas inlet temperature
        T_wi = inputs['T_cold_in']  # Water inlet temperature
        R1 = inputs['R_f_hot']  # Flue gas side fouling resistance
        R2 = inputs['R_f_cold']  # Water side fouling resistance
        
        # Step 2: Initialize guess
        T_fo_guess = T_fi - 10.0  # Initial guess for flue gas outlet (T_fo,guess <= T_fi)
        T_fo = T_fi  # Initialize T_fo = T_fi
        
        iteration = 0
        error = float('inf')
        
        # Step 3: Iterative solution
        while error > TOLERANCE and iteration < MAX_ITER:
            # (a) Evaluate flue gas properties at bulk mean temperature
            T_flue_mean = (T_fi + T_fo_guess) / 2
            props_flue = self.fluid_hot.get_properties(T_flue_mean)
            
            # (b) Compute flue-side fouling variables (Eqns 10-13)
            res_flue = self._calc_h_flue(m_dot_hot, props_flue, R1)
            h_flue = res_flue['h']
            
            # (c) Evaluate T_wo,guess from energy balance and water properties
            # Energy balance: Q = m_hot * cp_hot * (T_fi - T_fo_guess)
            Q_guess = m_dot_hot * props_flue['cp'] * (T_fi - T_fo_guess)
            
            # Estimate water outlet using inlet properties
            props_water_in = self.fluid_cold.get_properties(T_wi)
            T_wo_guess = T_wi + Q_guess / (m_dot_cold * props_water_in['cp'])
            
            # Water properties at bulk mean temperature
            T_water_mean = (T_wi + T_wo_guess) / 2
            props_water = self.fluid_cold.get_properties(T_water_mean)
            
            # (d) Determine water-side fouling variables (Eqns 4-7)
            res_water = self._calc_h_water(m_dot_cold, props_water, R2)
            h_water = res_water['h']
            
            # (e) Calculate U, NTU, epsilon
            U = self.calc_overall_U(h_flue, h_water, R1, R2)
            NTU = self.calc_ntu(U, m_dot_hot, props_flue['cp'], 
                                m_dot_cold, props_water['cp'])
            epsilon = self.calc_effectiveness(U, NTU, m_dot_hot, props_flue['cp'],
                                             m_dot_cold, props_water['cp'])
            
            # (f) Obtain T_fo and T_wo from computations
            C_h = m_dot_hot * props_flue['cp']
            C_c = m_dot_cold * props_water['cp']
            C_min = min(C_h, C_c)
            
            Q_actual = epsilon * C_min * (T_fi - T_wi)
            T_fo = T_fi - (Q_actual / C_h)
            T_wo = T_wi + (Q_actual / C_c)
            
            # (g) Update guesses using relaxation
            T_fo_guess = (T_fo_guess + T_fo) / 2
            T_wo_guess = (T_wo_guess + T_wo) / 2
            
            # Check convergence
            error = abs(T_fo_guess - T_fo)
            iteration += 1
        
        # Step 4: Record data
        return {
            'inputs': inputs,
            'results': {
                'T_hot_out': T_fo,  # T_fo
                'T_cold_out': T_wo,  # T_wo
                'Q_actual': Q_actual,
                'U_overall': U,
                'NTU': NTU,
                'epsilon': epsilon,
                'h_flue': h_flue,
                'h_water': h_water,
                'D_f_o': res_flue['D_f_o'],
                'D_f_i': res_water['D_f_i'],
                'V_max_f': res_flue['V_max_f'],
                'V_f': res_water['V_f'],
                'M_dot_flue_f': res_flue['M_dot_f'],
                'M_dot_water_f': res_water['M_dot_f'],
                'error': error,
                'iteration': iteration,
            }
        }

