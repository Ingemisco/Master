import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict

def parse_file(file_path: str) -> List[Tuple[float, int, List[int]]]:
    """
    Parse a single file and return the data as list of tuples:
    (float_value, list_length, list_content)
    """
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
                
            float_val = float(parts[0])
            n = int(parts[1])
            numbers = list(map(int, parts[2:2+n]))
            
            data.append((float_val, n, numbers))
    
    return data

def process_single_file(data: List[Tuple[float, int, List[int]]]) -> List[Tuple[float, float, int]]:
    """
    Process data from a single file to create intervals and normalize
    Returns: list of (start, end, length) tuples
    """
    if not data:
        return []
    
    # Sort by float value (descending order as in file)
    # data.sort(key=lambda x: x[0], reverse=True)
    
    # Get the maximum value for normalization
    max_val = data[0][0]
    
    intervals = []
    
    # First interval: [0, normalized_first_value]
    # first_val = data[0][0] / max_val
    # intervals.append((0.0, first_val, data[0][1]))
    
    # Middle intervals
    for i in range(1, len(data)):
        start = data[i][0] / max_val
        end = data[i-1][0] / max_val
        intervals.append((start, end, data[i][1]))
    
    # Last interval: [normalized_last_value, 1.0]
    # last_val = data[-1][0] / max_val
    # intervals.append((last_val, 1.0, data[-1][1]))
    
    return intervals

def process_directory(directory_path: str) -> Dict[float, List[int]]:
    """
    Process all files in a directory and collect length values for each normalized x position
    """
    # Collect all interval data from all files
    all_intervals = []
    
    for filename in os.listdir(directory_path):
        if filename.endswith('.txt') or not filename.startswith('.'):
            file_path = os.path.join(directory_path, filename)
            try:
                data = parse_file(file_path)
                intervals = process_single_file(data)
                all_intervals.append(intervals)
            except Exception as e:
                print(f"Error processing file {filename}: {e}")
    
    if not all_intervals:
        return {}
    
    # Create bins for averaging
    bin_edges = np.linspace(0, 1, 100)  # 100 bins from 0 to 1
    bin_lengths = {edge: [] for edge in bin_edges}
    
    # For each interval set, assign length values to bins
    for intervals in all_intervals:
        for start, end, length in intervals:
            # Find all bins that fall within this interval
            for bin_edge in bin_edges:
                if start <= bin_edge <= end:
                    bin_lengths[bin_edge].append(length)
    
    # Calculate averages
    averaged_data = {}
    for bin_edge, lengths in bin_lengths.items():
        if lengths:
            averaged_data[bin_edge] = np.mean(lengths)
        else:
            averaged_data[bin_edge] = 0
    
    return averaged_data

def plot_data(intervals: List[Tuple[float, float, int]], title: str = "List Lengths vs Normalized Intervals"):
    """
    Plot the step function for the given intervals
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for start, end, length in intervals:
        ax.hlines(y=length, xmin=start, xmax=end, linewidth=2, color='blue')
        # Add small vertical lines at interval boundaries for clarity
        ax.vlines(x=start, ymin=length-0.1, ymax=length+0.1, color='blue', linewidth=1, alpha=0.7)
        ax.vlines(x=end, ymin=length-0.1, ymax=length+0.1, color='blue', linewidth=1, alpha=0.7)
    
    ax.set_xlabel('Normalized Interval Value')
    ax.set_ylabel('List Length')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def plot_averaged_data(averaged_data: Dict[float, float], title: str = "Averaged List Lengths vs Normalized Intervals"):
    """
    Plot the averaged data as a step function
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Sort the data by x values
    x_values = sorted(averaged_data.keys())
    y_values = [averaged_data[x] for x in x_values]
    
    # Create step plot
    ax.step(x_values, y_values, where='post', linewidth=2, color='red')
    
    ax.set_xlabel('Normalized Interval Value')
    ax.set_ylabel('Average List Length')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot list lengths from data files')
    parser.add_argument('input_path', help='Path to input file or directory')
    args = parser.parse_args()
    
    input_path = args.input_path
    
    if os.path.isfile(input_path):
        # Process single file
        print(f"Processing file: {input_path}")
        data = parse_file(input_path)
        intervals = process_single_file(data)
        plot_data(intervals, f"List Lengths - {os.path.basename(input_path)}")
        
    elif os.path.isdir(input_path):
        # Process directory
        print(f"Processing directory: {input_path}")
        averaged_data = process_directory(input_path)
        plot_averaged_data(averaged_data, f"Averaged List Lengths - {os.path.basename(input_path)}")
        
    else:
        print(f"Error: {input_path} is not a valid file or directory")

if __name__ == "__main__":
    main()
