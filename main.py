"""
Main Simulation Script

Configurable heat exchanger simulation with support for:
- Multiple exchanger types (PHE, CrossFlow)
- Time column generation
- Noise addition
"""
import random
import pandas as pd
from datetime import timedelta, datetime
from config import (
    RANGES, SIMULATION_START_TIME, SIMULATION_END_TIME, TIME_FORMAT,
    EXCHANGER_TYPE, GENERATE_TIME_COLUMN, ADD_NOISE, NOISE_LEVEL, PERCENTAGE_NOISE, NOISE_COLUMNS,
    NUM_RUNS, HOURS_PER_RUN, FOULING_GROWTH_HOT, FOULING_GROWTH_COLD
)
from solver import Solver
from models.plate import PlateHeatExchanger
from models.crossflow import CrossFlowHeatExchanger
from utils.noise import add_noise_to_dataframe


def create_exchanger(exchanger_type: str = None):
    """
    Create heat exchanger instance based on type.
    
    Args:
        exchanger_type: 'PHE' or 'CROSSFLOW' (defaults to config value)
        
    Returns:
        Heat exchanger instance
    """
    if exchanger_type is None:
        exchanger_type = EXCHANGER_TYPE
    
    exchanger_type = exchanger_type.upper()
    
    if exchanger_type == 'PHE':
        return PlateHeatExchanger()
    elif exchanger_type == 'CROSSFLOW':
        return CrossFlowHeatExchanger()
    else:
        raise ValueError(f"Unknown exchanger type: {exchanger_type}. Must be 'PHE' or 'CROSSFLOW'")


def generate_random_inputs():
    """
    Step 1 from Table 2: Random Selection from Uniform Distribution
    
    Returns:
        Dictionary of random inputs within configured ranges
    """
    inputs = {}
    for key, (min_val, max_val) in RANGES.items():
        inputs[key] = random.uniform(min_val, max_val)
    return inputs


def generate_random_timestamp(start_date, end_date):
    """
    Generates a random timestamp between start_date and end_date
    
    Args:
        start_date: Start datetime
        end_date: End datetime
        
    Returns:
        Random timestamp between start and end
    """
    delta = end_date - start_date
    random_seconds = random.randint(0, int(delta.total_seconds()))
    return start_date + timedelta(seconds=random_seconds)


def generate_sequential_run(solver: Solver, start_time, hours: int = None,
                           fouling_growth_hot: tuple = None,
                           fouling_growth_cold: tuple = None,
                           run_id: int = None, end_time: datetime = None):
    """
    Simulates a continuous operation run with fouling growth.
    Fouling grows linearly or asymptotically with time.
    
    Args:
        solver: Solver instance
        start_time: Starting datetime for the run
        hours: Number of hours to simulate (defaults to config)
        fouling_growth_hot: Tuple (min, max) for hot side fouling growth per hour
        fouling_growth_cold: Tuple (min, max) for cold side fouling growth per hour
        run_id: Optional run identifier to track which run produced each row
        end_time: Optional end datetime - run will stop if current_time exceeds this
        
    Returns:
        List of data dictionaries
    """
    if hours is None:
        hours = HOURS_PER_RUN
    if fouling_growth_hot is None:
        fouling_growth_hot = FOULING_GROWTH_HOT
    if fouling_growth_cold is None:
        fouling_growth_cold = FOULING_GROWTH_COLD
    
    data_sequence = []
    current_time = start_time
    
    # Initial clean conditions
    current_Rf_hot = 0.0
    current_Rf_cold = 0.0
    
    for i in range(hours):
        # Advance time
        current_time += timedelta(hours=1)
        
        # Stop if we've exceeded the end time
        if end_time is not None and current_time > end_time:
            break
        
        # Fouling growth (simple linear growth model)
        current_Rf_hot += random.uniform(*fouling_growth_hot)
        current_Rf_cold += random.uniform(*fouling_growth_cold)
        
        # Operational fluctuations (sensor noise)
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
                **res['inputs'],
                **res['results']
            }
            
            # Add timestamp if configured
            if GENERATE_TIME_COLUMN:
                flat_data['timestamp'] = current_time.strftime(TIME_FORMAT)
            
            # Add run_id if provided
            if run_id is not None:
                flat_data['run_id'] = run_id
            
            data_sequence.append(flat_data)
            
        except Exception as e:
            print(f"Simulation failed at hour {i}: {e}")
            continue
    
    return data_sequence


def generate_random_dataset(solver: Solver, n_samples: int, with_time: bool = None):
    """
    Generate random dataset without sequential time dependency.
    
    Args:
        solver: Solver instance
        n_samples: Number of samples to generate
        with_time: Whether to include timestamp (defaults to config)
        
    Returns:
        List of data dictionaries
    """
    if with_time is None:
        with_time = GENERATE_TIME_COLUMN
    
    dataset = []
    
    for i in range(n_samples):
        # Random input generation
        inputs = generate_random_inputs()
        
        try:
            res = solver.run_simulation(inputs)
            
            flat_data = {
                **res['inputs'],
                **res['results']
            }
            
            # Add random timestamp if configured
            if with_time:
                timestamp = generate_random_timestamp(SIMULATION_START_TIME, SIMULATION_END_TIME)
                flat_data['timestamp'] = timestamp.strftime(TIME_FORMAT)
            
            dataset.append(flat_data)
            
        except Exception as e:
            print(f"Convergence failed for sample {i}: {e}")
            continue
    
    return dataset


def main(exchanger_type: str = None, with_time: bool = None, 
         with_noise: bool = None, n_samples: int = None, 
         sequential: bool = True, output_file: str = None):
    """
    Main simulation function.
    
    Args:
        exchanger_type: 'PHE' or 'CROSSFLOW' (defaults to config)
        with_time: Whether to include timestamp column (defaults to config)
        with_noise: Whether to add noise to output (defaults to config)
        n_samples: Number of samples for random generation (if sequential=False)
        sequential: Whether to use sequential fouling growth (default: True)
        output_file: Output CSV filename (default: auto-generated)
    """
    # Use config defaults if not provided
    if exchanger_type is None:
        exchanger_type = EXCHANGER_TYPE
    if with_time is None:
        with_time = GENERATE_TIME_COLUMN
    if with_noise is None:
        with_noise = ADD_NOISE
    
    # Create exchanger and solver
    exchanger = create_exchanger(exchanger_type)
    solver = Solver(exchanger)
    
    print(f"Generating {exchanger_type} heat exchanger simulation data...")
    print(f"  Time column: {with_time}")
    print(f"  Noise: {with_noise}")
    print(f"  Sequential: {sequential}")
    
    all_data = []
    
    if sequential:
        # Generate sequential fouling data
        total_hours_needed = NUM_RUNS * HOURS_PER_RUN
        available_hours = int((SIMULATION_END_TIME - SIMULATION_START_TIME).total_seconds() / 3600)
        
        print(f"Generating {NUM_RUNS} runs of {HOURS_PER_RUN} hours each...")
        print(f"  Total hours needed: {total_hours_needed}")
        print(f"  Available hours in date range: {available_hours}")
        
        if total_hours_needed > available_hours:
            print(f"  WARNING: Requested {total_hours_needed} hours but only {available_hours} hours available.")
            print(f"  Runs will be truncated at {SIMULATION_END_TIME.strftime(TIME_FORMAT)}")
        
        for run_idx in range(NUM_RUNS):
            # Shift start date for each run by HOURS_PER_RUN to prevent overlap
            # Each run lasts HOURS_PER_RUN hours, so space them accordingly
            run_start = SIMULATION_START_TIME + timedelta(hours=run_idx * HOURS_PER_RUN)
            
            # Skip runs that start after the end time
            if run_start >= SIMULATION_END_TIME:
                print(f"  Run {run_idx + 1}/{NUM_RUNS}: Skipped (starts after end time)")
                continue
            
            # Generate sequential run with run_id and end_time constraint
            run_data = generate_sequential_run(solver, run_start, run_id=run_idx, end_time=SIMULATION_END_TIME)
            all_data.extend(run_data)
            
            actual_hours = len(run_data)
            expected_hours = HOURS_PER_RUN
            if actual_hours < expected_hours:
                print(f"  Run {run_idx + 1}/{NUM_RUNS}: {actual_hours} samples (truncated from {expected_hours}, start: {run_start.strftime(TIME_FORMAT)})")
            else:
                print(f"  Run {run_idx + 1}/{NUM_RUNS}: {actual_hours} samples (start: {run_start.strftime(TIME_FORMAT)})")
    else:
        # Generate random dataset
        if n_samples is None:
            n_samples = NUM_RUNS * HOURS_PER_RUN
        
        print(f"Generating {n_samples} random samples...")
        all_data = generate_random_dataset(solver, n_samples, with_time)
    
    # Convert to DataFrame
    df = pd.DataFrame(all_data)
    
    # Add noise if configured
    if with_noise:
        print(f"Adding noise ({PERCENTAGE_NOISE*100:.1f}% of rows, level: {NOISE_LEVEL*100:.1f}% std dev)...")
        df = add_noise_to_dataframe(df, noise_level=NOISE_LEVEL, percentage_noise=PERCENTAGE_NOISE, noise_columns=NOISE_COLUMNS)
    
    # Generate output filename if not provided
    if output_file is None:
        exchanger_suffix = exchanger_type.lower()
        time_suffix = "_timestamp_" + datetime.now().strftime("%Y%m%d_%H%M%S") if with_time else ""
        noise_suffix = "_noisy" if with_noise else ""
        output_file = f"data/{exchanger_suffix}_simulation_dataset{time_suffix}{noise_suffix}.csv"
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    
    print(f"\nDataset generated: {output_file} ({len(df)} samples)")
    
    # Display sample
    display_cols = ['R_f_hot', 'T_hot_out', 'U_overall']
    if with_time:
        display_cols.insert(0, 'timestamp')
    print(f"\nSample data:")
    print(df[display_cols].head())
    
    return df


if __name__ == "__main__":
    # Run with config defaults
    main()
    
    # Example: Run with custom parameters
    # main(exchanger_type='CROSSFLOW', with_time=True, with_noise=True, sequential=False, n_samples=1000)
