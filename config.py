#Simulation Constraints based on literature
"""
Configuration module for heat exchanger fouling simulation.
This module defines all simulation constraints, physical constants, solver settings,
and noise/output configurations for modeling fouling behavior in plate heat exchangers (PHE).
Constants:
    RANGES (dict): Operating parameter ranges for simulation:
        - m_dot_hot (tuple): Hot side (flue gas) mass flow rate in kg/s
        - m_dot_cold (tuple): Cold side (water) mass flow rate in kg/s
        - T_hot_in (tuple): Hot side inlet temperature in Celsius
        - T_cold_in (tuple): Cold side inlet temperature in Celsius
        - R_f_hot (tuple): Hot side fouling resistance in m²·K/W
        - R_f_cold (tuple): Cold side fouling resistance in m²·K/W
    K_DEPOSIT (float): Thermal conductivity of fouling deposit in W/m·K
    K_PLATE (float): Thermal conductivity of plate material in W/m·K
    TOLERANCE (float): Convergence tolerance for iterative solver
    MAX_ITER (int): Maximum number of iterations for solver
    EXCHANGER_TYPE (str): Type of heat exchanger ('PHE' or 'CROSSFLOW')
    ADD_NOISE (bool): Flag to enable noise addition to output data
    PERCENTAGE_NOISE (float): Percentage noise level
    BASE_NOISE_LEVEL (float): Base noise level as fraction of measurement
    NOISE_LEVEL (float): Standard deviation of noise as fraction
    NOISE_COLUMNS (list): Column names to apply noise to
    RANDOM_SEED (int): Seed for reproducible random number generation
    NUM_RUNS (int): Number of separate maintenance cycles to simulate
    HOURS_PER_RUN (int): Hours of simulation data per run
    GENERATE_TIME_COLUMN (bool): Flag to generate timestamp column
    SIMULATION_START_DATE (datetime): Start date for simulation
    SIMULATION_START_TIME (datetime): Start datetime for simulation
    SIMULATION_END_TIME (datetime): End datetime for simulation
    TIME_FORMAT (str): Format string for timestamp output
    FOULING_GROWTH_HOT (tuple): Range for hot side fouling growth rate per hour
    FOULING_GROWTH_COLD (tuple): Range for cold side fouling growth rate per hour
    TAU_HOT (int): Time constant for hot side fouling growth in hours
    TAU_COLD (int): Time constant for cold side fouling growth in hours
"""
# Adjusted for the specific PHE geometry and fluid properties
from datetime import datetime

RANGES = {
    'm_dot_hot': (4.0, 8.5),        #kg/s --Unfouled flue gas flow rate ( ˙ M flue ) 
    'm_dot_cold': (0.75, 10.65),    #kg/s --Unfouled water flow rate ( ˙ M w )
    'T_hot_in': (170.0, 250),        #Celsius --flue inlet temperature ( T fi) 
    'T_cold_in': (20, 30),        #Celsius --water inlet temperature ( T wi 
    #fouling factor(m^2.K/W)
    'R_f_hot': (0.0, 0.00175),    #flue gas side fouling resistance ( R f,flue )
    'R_f_cold': (0.0, 0.00053),     #water side fouling resistance ( R f,water )
}


#Physical Constants
K_DEPOSIT = 2.3 #W/m.K
K_PLATE = 16.0 #W/m.K


#Solver Settings
TOLERANCE = 1e-6 #Convergence criteria
MAX_ITER = 100 #Maximum number of iterations

#Simulation Configuration
EXCHANGER_TYPE = 'CROSSFLOW'  # Options: 'PHE' or 'CROSSFLOW'




#NOISE CONFIGURATION
 # Whether to include timestamp column
ADD_NOISE = True  # Whether to add noise to output data
PERCENTAGE_NOISE = 0.15 # Noise level as percentage of measurement (e.g., 0.02 = 2%)
BASE_NOISE_LEVEL = 0.00  # Base noise level (0.02 = 2% standard deviation)
NOISE_LEVEL = 0.02  # Noise level as fraction (0.02 = 2% standard deviation)
NOISE_COLUMNS = [  # Columns to add noise to (if ADD_NOISE is True)
    'T_hot_out',
    'T_cold_out',
    'Q_actual',
    'U_overall',
    'h_hot',
    'h_cold',
]
RANDOM_SEED = 42  # For reproducibility of noise and random inputs
#Sequential Run Configuration
NUM_RUNS = 10  # Number of separate maintenance cycles/runs
HOURS_PER_RUN = 1000  # Hours of data per run


#Time Configurations
GENERATE_TIME_COLUMN = False
SIMULATION_START_DATE = datetime(2023, 1, 1)
SIMULATION_START_TIME = datetime(2023,1,1)
SIMULATION_END_TIME = datetime(2024,1,1)
TIME_FORMAT = "%Y-%m-%d %H:%M:%S"


#FOULING GROWTH CONFIGURATION
FOULING_GROWTH_HOT = (0.0, 0.00001)  # Range for fouling growth per hour (hot side)
FOULING_GROWTH_COLD = (0.0, 0.000005)  # Range for fouling growth per hour (cold side)
TAU_HOT = 1200  # Time constant for hot side fouling growth (hours)
TAU_COLD = 700  # Time constant for cold side fouling growth (hours)