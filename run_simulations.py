import os
import subprocess
import psutil
import time
import argparse
from tqdm import tqdm
from datetime import datetime

def get_quarter_from_month(month):
    """Convert month (1-12) to quarter (1-4)."""
    return (month - 1) // 3 + 1

def parse_args():
    parser = argparse.ArgumentParser(description='Run MCP simulations')
    parser.add_argument('--demand-type', choices=['gravity', 'db1b'], required=True,
                       help='Type of demand model to use')
    
    # Common parameters
    parser.add_argument('--year', type=int, required=True,
                       help='Year for the simulation')
    parser.add_argument('--month', type=int, required=True,
                       help='Month for flight schedule')
    parser.add_argument('--date', type=int, required=True,
                       help='Date for flight schedule')
    parser.add_argument('--coop-type', required=True,
                       help='Type of cooperation (e.g., alliance)')
    parser.add_argument('--cost-type', type=int, choices=[1, 2, 3], default=2,
                       help='Cost type: 1=distance, 2=time, 3=seats')
    
    # Gravity model specific parameters
    gravity_group = parser.add_argument_group('gravity model parameters')
    gravity_group.add_argument('--param1', type=float,
                               help='Population parameter 1 for gravity model')
    gravity_group.add_argument('--param2', type=float,
                               help='Population parameter 2 for gravity model')
    gravity_group.add_argument('--param3', type=float,
                               help='Decay parameter for gravity model')
    gravity_group.add_argument('--param4', type=int,
                               help='Distance parameter for gravity model')
    
    # DB1B specific parameters
    # db1b_group = parser.add_argument_group('db1b parameters')
    
    return parser.parse_args()

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
        demand_file = (f'{base_paths["derived"]}/gravity_demand_pop1_{kwargs["param1"]:.1f}_'
                      f'pop2_{kwargs["param2"]:.1f}_decay_{kwargs["param3"]:.1f}_d_{kwargs["param4"]}.dat')
        output_file = (f'{base_paths["results"]}/Y{kwargs["year"]}M{kwargs["month"]}D{kwargs["date"]}/results_gravity_pop1_{kwargs["param1"]:.1f}_'
                      f'pop2_{kwargs["param2"]:.1f}_decay_{kwargs["param3"]:.1f}_d_{kwargs["param4"]}_'
                      f'{kwargs["coop_type"]}_{cost_dict[cost_type]}.dat')
    
    elif demand_type == 'db1b':
        demand_file = f'{base_paths["derived"]}/DB1B_demand_Y{kwargs["year"]}_Q{get_quarter_from_month(kwargs["month"])}.dat'
        output_file = f'{base_paths["results"]}/Y{kwargs["year"]}M{kwargs["month"]}D{kwargs["date"]}/results_db1b_{kwargs["coop_type"]}_{cost_dict[cost_type]}.dat'
    
    else:
        raise ValueError(f"Unknown demand type: {demand_type}")
    
    return {
        'demand_file': demand_file,
        'flight_file': flight_file,
        'network_file': network_file,
        'output_file': output_file,
        'cost_type': cost_type
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
            # capture_output=True,
            # text=True
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



# def write_log(log_file, result):
#     timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
#     status = result['status']
#     duration = result['duration']
#     demand = result['demand']
#     flight = result['flight']
#     network = result['network']
#     cost_type = result['cost_type']
    
#     log_entry = f"{timestamp} | {demand} | {flight} | {network} | {cost_type} | {duration:8.2f}s | {status}"
    
#     if status == 'error':
#         log_entry += f" | Error: {result['error_msg']}"
    
#     with open(log_file, 'a') as f:
#         f.write(log_entry + '\n')
        
def main():
    args = parse_args()
    # Define directories
    input_dir = './data/derived'
    output_dir = './data/results'
    log_dir = './logs'
    
    # Create output and log directories if they don't exist
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    # Create log file with timestamp
    log_file = os.path.join(log_dir, f'simulation_log.txt')
    
    # Write header to log file
    with open(log_file, 'a') as f:
        f.write("timestamp | demand | flight | network | cost_type | duration | status\n")
        f.write("-" * 100 + "\n")

    if args.demand_type == 'gravity':
        if not all([args.param1, args.param2, args.param3, args.param4]):
            raise ValueError("All gravity model parameters are required")
        config = create_simulation_config(**vars(args))
    else:
        config = create_simulation_config(**vars(args))
    
    results = run_simulation(config)
    # write_log(log_file, results)
    

if __name__ == "__main__":
    main()