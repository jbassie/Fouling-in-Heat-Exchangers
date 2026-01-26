import random
import pandas as pd
from config import RANGES
from solver import Solver


def generate_randon_inputs():
    """
    Step 1 from Table 2: Random Selection from Uniform Distribution
    """
    inputs = {}
    for key, (min_val, max_val) in RANGES.items():
        inputs[key] = random.uniform(min_val, max_val)
    return inputs


def main():
    N_SAMPLES =  10000  #Number of Simulations
    solver = Solver()
    dataset = []

    print(f"Generating {N_SAMPLES} random input sets...")

    for i in range(N_SAMPLES):
        #1. Random Input Generation
        inp = generate_randon_inputs()

        #2. Run Iterative Solver
        try: 
            data_point = solver.run_simulation(inp)

            #Flatten structure for CSV
            flat_data = {**data_point['inputs'], **data_point['results']}
            dataset.append(flat_data)
        except Exception as e:
            print(f"Covergence failed for sample {i}: {e}")
    
    df = pd.DataFrame(dataset)
    df.to_csv('phe_simulation_dataset.csv', index=False)
    print(f"Dataset saved to phe_simulation_dataset.csv with {len(dataset)} samples")
    print(df.head())


if __name__ == '__main__':
    main()
