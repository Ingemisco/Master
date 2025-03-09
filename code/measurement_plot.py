import os
import json
import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.stats import linregress

def regression(x, y):
    log_x = np.log(x)
    log_y = np.log(y)

    slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)

    k = slope  # Exponent
    a = np.exp(intercept)  # Coefficient

    return a, k

def read_json_files(directory):
    data_points = []
    
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".json"):
            with open(os.path.join(directory, filename), 'r') as f:
                data = json.load(f)
                data_points.append(data)
    
    return data_points

def plot_data(directory):
    data_points = read_json_files(directory)
    
    if not data_points:
        print("No JSON files found!")
        return
    
    # Extract values
    min_times = [entry["min"] for entry in data_points]
    max_times = [entry["max"] for entry in data_points]
    avg_times = [entry["avg"] for entry in data_points]
    
    x_labels =  [entry["min_point_count"] for entry in data_points] 

    _n = len(x_labels)
    coefficient, exponent = regression(x_labels[_n // 2:], avg_times[_n // 2:])
    regression_curve = [coefficient * (x ** exponent) for x in x_labels]
    
    # Extract metadata (same for all files)
    algorithm = data_points[0]["algorithm"]
    point_count = data_points[0]["min_point_count"]  # Same as max_point_count
    dimension_count = data_points[0]["min_dimension_count"]  # Same as max_dimension_count
    
    # Plot data
    plt.figure(figsize=(10, 5))
    plt.plot(x_labels, min_times, marker='o', label='Min Time')
    plt.plot(x_labels, max_times, marker='s', label='Max Time')
    plt.plot(x_labels, avg_times, marker='d', label='Avg Time')
    plt.plot(x_labels, regression_curve, marker='d', label=f"Curve y = {coefficient} * x^{exponent}")
    
    # Labels and title
    plt.xlabel("Polyline Length")
    plt.ylabel("Time (seconds)")
    plt.title(f"Performance of {algorithm}\nDimensions: {dimension_count}")
    plt.legend()
    plt.grid()
    
    # Show plot
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory>")
        sys.exit(1)
    
    plot_data(sys.argv[1])



