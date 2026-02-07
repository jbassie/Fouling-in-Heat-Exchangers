import random
import pandas as pd
from datetime import timedelta
from config import RANGES, SIMULATION_START_TIME, SIMULATION_END_TIME, TIME_FORMAT
from solver import Solver


def generate_randon_inputs():
    """
    Step 1 from Table 2: Random Selection from Uniform Distribution
    """
    inputs = {}
    for key, (min_val, max_val) in RANGES.items():
        inputs[key] = random.uniform(min_val, max_val)
    return inputs


def generate_random_timestamp(start_date, end_date):
    """
    Generates a random timestamp between start_date and end_date
    """
    delta = end_date - start_date
    random_seconds = random.randint(0, int(delta.total_seconds()))
    return start_date + timedelta(seconds=random_seconds)


def generate_sequential_run(solver, start_time, hours=24):
    """
    Alternative generation: Simulates a continuous operation run.
    Fouling grows linearly or asymptotically with time[cite: 119].
    """
    data_sequence = []
    current_time = start_time


    #Initial clean conditions
    current_Rf_hot = 0.0
    current_Rf_cold = 0.0

    for i in range(hours):
        #Generate random inputs for each hour
        current_time += timedelta(hours=1)

        #Fouling Grows(Simple linear growth model for demonstration)
        current_Rf_hot += random.uniform(0, 0.00001)  # Simulate fouling growth
        current_Rf_cold += random.uniform(0.0, 0.000005)  # Simulate fouling growth


        # Operational Fluctuations (Sensor Noise)
        # Inputs vary narrowly to simulate steady operation
        inputs = {
            'm_dot_hot': random.uniform(4.0, 4.5),
            'm_dot_cold': random.uniform(2.0, 2.5),
            'T_hot_in': random.uniform(195, 205),
            'T_cold_in': random.uniform(20, 25),
            'R_f_hot': current_Rf_hot,
            'R_f_cold': current_Rf_cold,
        }
        try:
            res = solver.run_simulation(inputs)

            flat_data = {
                'timestamp': current_time.strftime(TIME_FORMAT),
                **res['inputs'],
                **res['results']
            }
            data_sequence.append(flat_data)
            
        except Exception as e:
            continue
        # Record
        
       
    return data_sequence


def main():
      #Number of Simulations
    solver = Solver()
    all_data = []

    print("Generating sequential fouling data for Transformer training...")
    
    # Generate 10 separate 'runs' (e.g., 10 different maintenance cycles)
    for run_idx in range(10):
        # Shift start date for each run
        run_start = SIMULATION_START_TIME + timedelta(days=run_idx*30)
        
        # Generate 500 hours of data per run
        run_data = generate_sequential_run(solver, run_start, hours=500)
        all_data.extend(run_data)

    # Save to CSV
    df = pd.DataFrame(all_data)
    df.to_csv('phe_transformer_dataset_timestamp.csv', index=False)
    
    print(f"Dataset generated: phe_transformer_dataset_timestamp.csv ({len(df)} samples)")
    print(df[['timestamp', 'R_f_hot', 'T_hot_out', 'U_overall']].head())

if __name__ == "__main__":
    main()





# def main():
#     N_SAMPLES =  10000  #Number of Simulations
#     solver = Solver()
#     dataset = []

#     print(f"Generating {N_SAMPLES} random input sets...")

#     for i in range(N_SAMPLES):
#         #1. Random Input Generation
#         inp = generate_randon_inputs()

#         #2. Run Iterative Solver
#         try: 
#             data_point = solver.run_simulation(inp)

#             #Flatten structure for CSV
#             flat_data = {**data_point['inputs'], **data_point['results']}
#             dataset.append(flat_data)
#         except Exception as e:
#             print(f"Covergence failed for sample {i}: {e}")
    
#     df = pd.DataFrame(dataset)
#     df.to_csv('phe_simulation_dataset.csv', index=False)
#     print(f"Dataset saved to phe_simulation_dataset.csv with {len(dataset)} samples")
#     print(df.head())
