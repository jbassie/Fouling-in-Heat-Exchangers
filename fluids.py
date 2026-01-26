"""
Handles fluid property retrival.(Simplified linear approximations for demonstration purposes)

"""

class FluidProperties:
    """
    Evaluate fluid properties at bulk mean temperature.
    Logic corresponds to step 3a and 3c in literature.
    """

    def __init__(self, fluid_type):
        self.fluid_type = fluid_type

    def get_properties(self, T_Celcuis):
        """
        Return dictionary of fluid properties at given temperature.
        {
        'rho': density (kg/m^3),
        'mu': viscosity (Pa.s),
        'cp': specific heat (J/kg.K),
        'k': thermal conductivity (W/m.K),
        'pr': Prandtl number,
        }

        """
        #Simplified properties for Air (Hot Fluid) and Water (Cold Fluid)
        # In a real scenerio, use CoolProp or Steam Tables
        if self.fluid_type == 'flue_gas':
            #Air-like properties at 200C
            return {
                'rho': 0.746,
                'mu': 2.5e-5,
                'cp': 1026.0,
                'k': 0.0039,
                'pr': 0.68,
            }
        elif self.fluid_type == 'water':
            #Water-like properties at 30C
            return {
                'rho': 995.7,
                'mu': 0.0008,
                'cp': 4178.0,
                'k': 0.6,
                'pr': 5.42,
            }
        else:
            raise ValueError(f"Invalid fluid type: {self.fluid_type}")