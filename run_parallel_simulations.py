import os
import glob
import re
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from itertools import product
import time
from datetime import datetime

def get_quarter_from_month(month):
    """Convert month (1-12) to quarter (1-4)."""
    return (month - 1) // 3 + 1

def get_next_realization_number(base_output_path):
    """
    Get the next realization number by checking existing files.
    
    Args:
        base_output_path (str): Base path without realization number
        
    Returns:
        int: Next available realization number
    """
    # Look for existing files with pattern base_output_path_r{number}.dat
    existing_files = glob.glob(f"{base_output_path}_r[0-9]*.dat")
    if not existing_files:
        return 1
    
    # Extract realization numbers and find the maximum
    real_numbers = [int(re.search(r'_r(\d+)\.dat$', f).group(1)) for f in existing_files]
    return max(real_numbers) + 1

def create_simulation_config(demand_type='gravity', **kwargs):
    """
    Create simulation configuration based on demand type and parameters.
    
    Args:
        demand_type (str): Either 'gravity' or 'db1b'
        **kwargs: Parameters depending on demand_type
            For gravity:
                param1, param2, param3, param4: gravity model parameters
                year, month, date: date parameters
                coop_type: cooperation type
            For db1b:
                year, quarter: time parameters
                month, date: flight schedule date
                coop_type: cooperation type
    
    Returns:
        dict: Configuration with all file paths and parameters
    """
    base_paths = {
        'derived': './data/derived',
        'results': './data/results'
    }
    
    cost_type = kwargs.get('cost_type', 2)  # default to time-based cost
    cost_dict = {1: 'distance', 2: 'time', 3: 'seats'}
    
    # Common files
    flight_file = f'{base_paths["derived"]}/list_of_super_flights_Y{kwargs["year"]}_M{kwargs["month"]}_D{kwargs["date"]}.dat'
    network_file = f'{base_paths["derived"]}/fcn_Y{kwargs["year"]}_M{kwargs["month"]}_D{kwargs["date"]}_{kwargs["coop_type"]}.dat'
    
    if demand_type == 'gravity':
        demand_file = f'gravity_demand_pop1_{kwargs["param1"]:.1f}_pop2_{kwargs["param2"]:.1f}_decay_{kwargs["param3"]:.1f}_d_{kwargs["param4"]}.dat'
        demand_file = os.path.join(base_paths["derived"], demand_file)

        base_output = f'results_gravity_{kwargs["coop_type"]}_pop1_{kwargs["param1"]:.1f}_pop2_{kwargs["param2"]:.1f}_decay_{kwargs["param3"]:.1f}_d_{kwargs["param4"]}_{cost_dict[cost_type]}'
        base_output_path = os.path.join(base_paths["results"], f'Y{kwargs["year"]}M{kwargs["month"]}D{kwargs["date"]}', base_output)
    
    elif demand_type == 'db1b':
        demand_file = f'DB1B_demand_Y{kwargs["year"]}_Q{get_quarter_from_month(kwargs["month"])}.dat'
        demand_file = os.path.join(base_paths["derived"], demand_file)

        base_output = f'results_db1b_{kwargs["coop_type"]}_{cost_dict[cost_type]}'
        base_output_path = os.path.join(base_paths["results"], f'Y{kwargs["year"]}M{kwargs["month"]}D{kwargs["date"]}', base_output)
    
    else:
        raise ValueError(f"Unknown demand type: {demand_type}")
    
    # Create results directory if it doesn't exist
    results_dir = os.path.join(base_paths["results"], f'Y{kwargs["year"]}M{kwargs["month"]}D{kwargs["date"]}')
    os.makedirs(results_dir, exist_ok=True)

    # Get next realization number and create output file name
    realization = get_next_realization_number(base_output_path)
    output_file = f"{base_output_path}_r{realization}.dat"
    
    return {
        'demand_file': demand_file,
        'flight_file': flight_file,
        'network_file': network_file,
        'output_file': output_file,
        'cost_type': cost_type,
        'realization': realization
    }

def run_simulation(config):
    demand_file = config['demand_file']
    flight_file = config['flight_file']
    network_file = config['network_file']
    output_file = config['output_file']
    cost_type = config['cost_type']

    start_time = time.time()  # Track start time

    try:
        cmd = ['./min_cost_percolation.out', demand_file, flight_file, network_file, output_file, str(cost_type)]
        print(f"Executing: {' '.join(cmd)}")
        
        # Use shell=True to execute exactly as command line
        result = subprocess.run(
            ' '.join(cmd),  # Join command into a single string
            shell=True,     # Use shell to execute
            check=True,     # Raise exception on error
        )

        end_time = time.time()  # Track end time
        duration = end_time - start_time
        
        return {
            'demand': demand_file,
            'flight': flight_file,
            'network': network_file,
            'cost_type': cost_type,
            'duration': duration,
            'status': 'success'
        }
        
    except subprocess.CalledProcessError as e:
        end_time = time.time()
        duration = end_time - start_time
        print(f"Command failed with return code {e.returncode}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return {
            'demand': demand_file,
            'flight': flight_file,
            'network': network_file,
            'cost_type': cost_type,
            'duration': duration,
            'status': 'error',
            'error_msg': str(e)
        }


def parse_args():
    parser = argparse.ArgumentParser(description='Run multiple MCP simulations in parallel')
    parser.add_argument('--year', type=int, required=True,
                       help='Year for the simulation')
    parser.add_argument('--month', type=int, required=True,
                       help='Month for flight schedule')
    parser.add_argument('--date', type=int, required=True,
                       help='Date for flight schedule')
    parser.add_argument('--threads', type=int, default=3,
                       help='Number of parallel threads to run (default: 3)')
    parser.add_argument('--iterations', type=int, default=1,
                       help='Number of iterations to run for each parameter combination (default: 1)')    
    return parser.parse_args()


def get_simulation_parameters():
    """Define all parameter combinations to run"""
    # Define parameter ranges
    # demand_types = ['gravity', 'db1b']
    demand_types = ['gravity']
    coop_types = ['coop', 'non_coop']  # Add more cooperation types as needed
    cost_types = [1, 2, 3]
    
    # Gravity model specific parameters
    gravity_params = [
        # param1, param2, param3, param4
        (1.0, 1.0, 2, 300),
        # Add more parameter combinations as needed
    ]

    # Generate all combinations
    all_params = []
    
    # Add gravity model combinations
    for coop_type, cost_type, (param1, param2, param3, param4) in product(
        coop_types, cost_types, gravity_params):
        all_params.append({
            'demand_type': 'gravity',
            'coop_type': coop_type,
            'cost_type': cost_type,
            'param1': param1,
            'param2': param2,
            'param3': param3,
            'param4': param4
        })
    
    # Add DB1B combinations
    for coop_type, cost_type in product(coop_types, cost_types):
        all_params.append({
            'demand_type': 'db1b',
            'coop_type': coop_type,
            'cost_type': cost_type
        })
    
    return all_params

def write_log(log_file, params, result):
    """Write simulation results to log file"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    status = result['status']
    duration = result['duration']
    
    # Create detailed parameter string
    param_str = f"demand_type={params['demand_type']}"
    param_str += f" coop_type={params['coop_type']}"
    param_str += f" cost_type={params['cost_type']}"
    if params['demand_type'] == 'gravity':
        param_str += f" params={params['param1']},{params['param2']},{params['param3']},{params['param4']}"
    
    log_entry = f"{timestamp} | {param_str} | {duration:8.2f}s | {status}"
    
    if status == 'error':
        log_entry += f" | Error: {result['error_msg']}"
    
    with open(log_file, 'a') as f:
        f.write(log_entry + '\n')

def run_simulation_with_params(params, year, month, date, log_file):
    """Wrapper function to run a single simulation with given parameters"""
    start_time = time.time()
    
    try:
        # Combine temporal parameters with simulation parameters
        full_params = {
            **params,
            'year': year,
            'month': month,
            'date': date
        }
        
        # Create configuration and run simulation
        config = create_simulation_config(**full_params)
        simulation_result = run_simulation(config)
        
        result = {
            'status': 'success',
            'duration': time.time() - start_time,
            'data': simulation_result
        }
    except Exception as e:
        result = {
            'status': 'error',
            'duration': time.time() - start_time,
            'error_msg': str(e)
        }
    
    # Write to log
    write_log(log_file, params, result)
    return result

def main():
    args = parse_args()
    
    # Create log directory if it doesn't exist
    log_dir = './logs'
    os.makedirs(log_dir, exist_ok=True)
    
    # Create log file with timestamp
    log_file = os.path.join(log_dir, f'parallel_simulation_log_{args.year}_{args.month}_{args.date}.txt')
    
    # Write header to log file
    with open(log_file, 'w') as f:
        f.write("timestamp | parameters | duration | status\n")
        f.write("-" * 100 + "\n")
    
    # Get all parameter combinations
    all_params = get_simulation_parameters()
    all_params_with_iterations = all_params * args.iterations
    
    print(f"Starting {len(all_params)} simulations with {args.threads} parallel threads...")
    
    # Run simulations in parallel with user-specified number of threads
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(
                run_simulation_with_params,
                params,
                args.year,
                args.month,
                args.date,
                log_file
            )
            for params in all_params_with_iterations
        ]
        
        # Wait for all simulations to complete
        completed = 0
        for future in futures:
            result = future.result()
            completed += 1
            print(f"Progress: {completed}/{len(all_params_with_iterations)} simulations completed")
    
    print(f"\nAll simulations completed. Results logged to: {log_file}")

if __name__ == "__main__":
    main()
