#Simulation Constraints based on literature
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

#Time Configurations
SIMULATION_START_DATE = datetime(2023, 1, 1)
SIMULATION_START_TIME = datetime(2023,1,1)
SIMULATION_END_TIME = datetime(2024,1,1)
TIME_FORMAT = "%Y-%m-%d %H:%M:%S"

#Simulation Configuration
EXCHANGER_TYPE = 'PHE'  # Options: 'PHE' or 'CROSSFLOW'
GENERATE_TIME_COLUMN = True  # Whether to include timestamp column
ADD_NOISE = True  # Whether to add noise to output data
PERCENTAGE_NOISE = 0.15  # Percentage of noise to add to output data
NOISE_LEVEL = 0.02  # Noise level as fraction (0.02 = 2% standard deviation)
NOISE_COLUMNS = [  # Columns to add noise to (if ADD_NOISE is True)
    'T_hot_out',
    'T_cold_out',
    'Q_actual',
    'U_overall',
    'h_hot',
    'h_cold',
]

#Sequential Run Configuration
NUM_RUNS = 10  # Number of separate maintenance cycles/runs
HOURS_PER_RUN = 1000  # Hours of data per run
FOULING_GROWTH_HOT = (0.0, 0.00001)  # Range for fouling growth per hour (hot side)
FOULING_GROWTH_COLD = (0.0, 0.000005)  # Range for fouling growth per hour (cold side)