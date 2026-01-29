import random
import pandas as pd
from datetime import timedelta
from config import SIMULATION_START_DATE, TIME_FORMAT
from solver import Solver

def generate_sequential_run(solver, start_time, hours=24):
    """
    Simulates a continuous operation run.
    [cite_start]Fouling grows linearly or asymptotically with time[cite: 119].
    """
    data_sequence = []
    current_time = start_time

    # Initial clean conditions
    current_Rf_hot = 0.0
    current_Rf_cold = 0.0

    for i in range(hours):
        # 1. Advance Time
        current_time += timedelta(hours=1)

        # 2. Simulate Fouling Growth (Linear model)
        # Random increment mimics varying deposition rates
        current_Rf_hot += random.uniform(0, 0.00001) 
        current_Rf_cold += random.uniform(0, 0.000005)

        # 3. Operational Fluctuations (Sensor Noise)
        # Inputs vary narrowly to simulate steady operation
        inputs = {
            'm_dot_hot': random.uniform(4.0, 4.5),
            'm_dot_cold': random.uniform(2.0, 2.5),
            'T_hot_in': random.uniform(195, 205),
            'T_cold_in': random.uniform(20, 25),
            # CRITICAL FIX: Keys must match config.py and Solver expectations
            'R_f_hot': current_Rf_hot,   
            'R_f_cold': current_Rf_cold, 
        }

        try:
            # 4. Solve for this hour
            res = solver.run_simulation(inputs)

            # 5. Flatten Data
            flat_data = {
                'timestamp': current_time.strftime(TIME_FORMAT),
                **res['inputs'],
                **res['results']
            }
            data_sequence.append(flat_data)
            
        except Exception as e:
            print(f"Simulation failed at hour {i}: {e}")
            # Skip failed convergences (rare in steady operation)
            continue
       
    return data_sequence

def main():
    solver = Solver()
    all_data = []

    print("Generating sequential fouling data for Transformer training...")
    
    # Generate 10 separate 'runs' (e.g., 10 different maintenance cycles)
    for run_idx in range(10):
        # Shift start date for each run (e.g., every 30 days)
        run_start = SIMULATION_START_DATE + timedelta(days=run_idx*30)
        
        # Generate 500 hours of data per run
        run_data = generate_sequential_run(solver, run_start, hours=500)
        all_data.extend(run_data)

    # Save to CSV
    df = pd.DataFrame(all_data)
    filename = 'phe_transformer_dataset_timestamp.csv'
    df.to_csv(filename, index=False)
    
    print(f"Dataset generated: {filename} ({len(df)} samples)")
    
    # Validation Print
    # Checks specific columns to ensure keys were propagated correctly
    print(df[['timestamp', 'R_f_hot', 'T_hot_out', 'U_overall']].head())

if __name__ == "__main__":
    main()