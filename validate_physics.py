"""
Physics Validation Script for PHE Fouling Model

This script validates that the fouling model behaves correctly:
- Increasing fouling resistance should decrease U_overall and Q_actual
- Increasing fouling should cause T_hot_out to increase (hot side cools less)
- Increasing fouling should cause T_cold_out to decrease (cold side heats less)
"""
from solver import Solver


def validate_fouling_physics():
    """
    Test that increasing fouling resistance correctly reduces heat transfer.
    
    Returns:
        bool: True if all physics checks pass, False otherwise
    """
    solver = Solver()
    
    # Base inputs (held constant)
    base_inputs = {
        'm_dot_hot': 5.0,
        'm_dot_cold': 3.0,
        'T_hot_in': 200.0,
        'T_cold_in': 25.0,
    }
    
    # Test increasing R_f_hot (cold side clean)
    print("Testing increasing R_f_hot (R_f_cold = 0)...")
    results_hot = []
    
    for R_f_hot in [0.0, 0.0005, 0.001, 0.0015]:
        inputs = {**base_inputs, 'R_f_hot': R_f_hot, 'R_f_cold': 0.0}
        try:
            res = solver.run_simulation(inputs)
            results_hot.append({
                'R_f_hot': R_f_hot,
                'U_overall': res['results']['U_overall'],
                'Q_actual': res['results']['Q_actual'],
                'T_hot_out': res['results']['T_hot_out'],
                'T_cold_out': res['results']['T_cold_out'],
            })
        except Exception as e:
            print(f"  Error at R_f_hot={R_f_hot}: {e}")
            return False
    
    # Test increasing R_f_cold (hot side clean)
    print("Testing increasing R_f_cold (R_f_hot = 0)...")
    results_cold = []
    
    for R_f_cold in [0.0, 0.0001, 0.0002, 0.0003]:
        inputs = {**base_inputs, 'R_f_hot': 0.0, 'R_f_cold': R_f_cold}
        try:
            res = solver.run_simulation(inputs)
            results_cold.append({
                'R_f_cold': R_f_cold,
                'U_overall': res['results']['U_overall'],
                'Q_actual': res['results']['Q_actual'],
                'T_hot_out': res['results']['T_hot_out'],
                'T_cold_out': res['results']['T_cold_out'],
            })
        except Exception as e:
            print(f"  Error at R_f_cold={R_f_cold}: {e}")
            return False
    
    # Validate results
    print("\n=== Validation Results ===")
    
    # Check R_f_hot effects
    print("\nR_f_hot increasing:")
    all_pass = True
    
    for i in range(1, len(results_hot)):
        prev = results_hot[i-1]
        curr = results_hot[i]
        
        U_decreased = curr['U_overall'] < prev['U_overall']
        Q_decreased = curr['Q_actual'] < prev['Q_actual']
        T_hot_increased = curr['T_hot_out'] > prev['T_hot_out']
        T_cold_decreased = curr['T_cold_out'] < prev['T_cold_out']
        
        print(f"  R_f_hot: {prev['R_f_hot']:.6f} -> {curr['R_f_hot']:.6f}")
        print(f"    U: {prev['U_overall']:.2f} -> {curr['U_overall']:.2f} ({'PASS' if U_decreased else 'FAIL'})")
        print(f"    Q: {prev['Q_actual']:.1f} -> {curr['Q_actual']:.1f} ({'PASS' if Q_decreased else 'FAIL'})")
        print(f"    T_hot_out: {prev['T_hot_out']:.2f} -> {curr['T_hot_out']:.2f} ({'PASS' if T_hot_increased else 'FAIL'})")
        print(f"    T_cold_out: {prev['T_cold_out']:.2f} -> {curr['T_cold_out']:.2f} ({'PASS' if T_cold_decreased else 'FAIL'})")
        
        if not (U_decreased and Q_decreased and T_hot_increased and T_cold_decreased):
            all_pass = False
            print(f"    WARNING: Physics check FAILED for this step!")
    
    # Check R_f_cold effects
    print("\nR_f_cold increasing:")
    
    for i in range(1, len(results_cold)):
        prev = results_cold[i-1]
        curr = results_cold[i]
        
        U_decreased = curr['U_overall'] < prev['U_overall']
        Q_decreased = curr['Q_actual'] < prev['Q_actual']
        T_hot_increased = curr['T_hot_out'] > prev['T_hot_out']
        T_cold_decreased = curr['T_cold_out'] < prev['T_cold_out']
        
        print(f"  R_f_cold: {prev['R_f_cold']:.6f} -> {curr['R_f_cold']:.6f}")
        print(f"    U: {prev['U_overall']:.2f} -> {curr['U_overall']:.2f} ({'PASS' if U_decreased else 'FAIL'})")
        print(f"    Q: {prev['Q_actual']:.1f} -> {curr['Q_actual']:.1f} ({'PASS' if Q_decreased else 'FAIL'})")
        print(f"    T_hot_out: {prev['T_hot_out']:.2f} -> {curr['T_hot_out']:.2f} ({'PASS' if T_hot_increased else 'FAIL'})")
        print(f"    T_cold_out: {prev['T_cold_out']:.2f} -> {curr['T_cold_out']:.2f} ({'PASS' if T_cold_decreased else 'FAIL'})")
        
        if not (U_decreased and Q_decreased and T_hot_increased and T_cold_decreased):
            all_pass = False
            print(f"    WARNING: Physics check FAILED for this step!")
    
    print("\n" + "="*50)
    if all_pass:
        print("SUCCESS: All physics checks PASSED!")
        print("Fouling model behaves correctly.")
    else:
        print("ERROR: Some physics checks FAILED!")
        print("Fouling model needs further investigation.")
    
    return all_pass


if __name__ == "__main__":
    validate_fouling_physics()

