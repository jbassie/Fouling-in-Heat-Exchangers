#Simulation Constraints based on literature
# Adjusted for the specific PHE geometry and fluid properties

RANGES = {
    'm_dot_hot': (4.0, 8.5),        #kg/s
    'm_dot_cold': (0.75, 10.65),    #kg/s
    'T_hot_in': (170.0, 250),        #Celsius
    'T_cold_in': (20, 30),        #Celsius   #TODO: add a range for the cold fluid
    #fouling factor(m^2.K/W)
    'R_f_hot': (0.0, 0.00175),
    'R_f_cold': (0.0, 0.00053),
}

#Physical Constants
K_DEPOSIT = 2.3 #W/m.K
K_PLATE = 16.0 #W/m.K


#Solver Settings
TOLERANCE = 1e-6 #Convergence criteria
MAX_ITER = 100 #Maximum number of iterations