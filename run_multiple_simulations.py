import itertools
import subprocess
from pathlib import Path

def run_sweep():
    # Define parameter spaces
    params = {
        'demand_type': ['db1b', 'gravity'],
        'coop_type': ['coop', 'non_coop'],
        'cost_type': [1, 2, 3],
        'year': [2023],
        'month': [4],
        'date': [18]
    }
    
    # Gravity specific parameters
    gravity_params = {
        'param1': [1.0],
        'param2': [1.0],
        'param3': [2.0],
        'param4': [300]
    }

    # Generate all combinations
    base_keys = ['demand_type', 'coop_type', 'cost_type', 'year', 'month', 'date']
    base_values = [params[k] for k in base_keys]
    
    for values in itertools.product(*base_values):
        param_dict = dict(zip(base_keys, values))
        
        if param_dict['demand_type'] == 'gravity':
            # Add gravity parameter combinations
            gravity_keys = list(gravity_params.keys())
            gravity_values = [gravity_params[k] for k in gravity_keys]
            
            for g_values in itertools.product(*gravity_values):
                gravity_dict = dict(zip(gravity_keys, g_values))
                full_params = {**param_dict, **gravity_dict}
                run_simulation(full_params)
        else:
            run_simulation(param_dict)

def run_simulation(params):
    cmd = ['python3', 'run_simulations.py']
    for key, value in params.items():
        cmd.extend([f'--{key.replace("_", "-")}', str(value)])
    
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd)

if __name__ == "__main__":
    run_sweep()