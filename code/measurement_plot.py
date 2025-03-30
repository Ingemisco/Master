import os
import json
import matplotlib.pyplot as plt
from matplotlib.scale import FuncScale
import sys
import numpy as np
from scipy.stats import linregress

# Define forward and inverse cubic transformations
def cubic_transform(x):
    return np.sign(x) * np.abs(x)**3  # Preserves sign for negative values

def inverse_cubic_transform(x):
    return np.sign(x) * np.abs(x)**(1/3)  # Inverse cubic root

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
        if filename.endswith(".json") and filename != "specification.json":
            with open(os.path.join(directory, filename), 'r') as f:
                data = json.load(f)
                data_points.append(data)
    
    return data_points

def plot_data(directory):
    data_points = read_json_files(directory)
    
    if not data_points:
        print("No JSON files found!")
        return
    
    data_points.sort(key=lambda x: int(x["point_count"]))
    # Extract values
    min_times = [entry["min"] for entry in data_points]
    max_times = [entry["max"] for entry in data_points]
    avg_times = [entry["avg"] for entry in data_points]
    
    x_labels =  [entry["point_count"] for entry in data_points] 

    _n = len(x_labels)
    coefficient, exponent = regression(x_labels[_n // 4:], avg_times[_n // 4:])
    x = np.linspace(x_labels[0], x_labels[-1], 100)
    regression_curve = [coefficient * (_x ** exponent) for _x in x]
    
    # Extract metadata (same for all files)
    algorithm = data_points[0]["algorithm"]
    point_count = data_points[0]["point_count"]
    dimension_count = data_points[0]["dimension"]
    
    # Plot data
    plt.figure(figsize=(10, 5))
    plt.plot(x_labels, min_times, marker='o', label='Min Time')
    plt.plot(x_labels, max_times, marker='s', label='Max Time')
    plt.plot(x_labels, avg_times, marker='d', label='Avg Time')

    plt.plot(x, regression_curve, label=f"Curve y = {coefficient:.2e} * x^{exponent:.2f}")
    
    # Labels and title

    plt.xlabel("Polyline Length")
    plt.ylabel("Time (seconds)")
    plt.title(f"Performance of {algorithm}\nDimensions: {dimension_count}")
    plt.legend()
    ax = plt.gca()
    ax.grid(True)
    ax.set_yscale('function', functions=(cubic_transform, inverse_cubic_transform))
    time_ticks = np.arange(0, max_times[-1] + 0.21, 0.1)
    time_ticks = [inverse_cubic_transform(t) for t in time_ticks]
    ax.set_yticks(time_ticks)
    # ax.set_yscale('log')
    plt.grid()
    
    # Show plot
    plt.show()


def plot_data_compare(directory1, directory2):
    data_points1 = read_json_files(directory1)
    data_points2 = read_json_files(directory2)
    
    if not data_points1:
        print(f"No JSON files found in {directory1}!")
        return
    if not data_points2:
        print(f"No JSON files found in {directory2}!")
        return

    data_points1.sort(key=lambda x: int(x["point_count"]))
    data_points2.sort(key=lambda x: int(x["point_count"]))

    if [entry["point_count"] for entry in data_points1] != [entry["point_count"] for entry in data_points2]:
        print(f"Not same polyline sizes in the data")
        return

    # Extract values
    avg_times1 = [entry["avg"] for entry in data_points1]
    avg_times2 = [entry["avg"] for entry in data_points2]
    
    x_labels =  [entry["point_count"] for entry in data_points1]  

    _n = len(x_labels)
    coefficient1, exponent1 = regression(x_labels[_n // 4:], avg_times1[_n // 4:])
    coefficient2, exponent2 = regression(x_labels[_n // 4:], avg_times2[_n // 4:])
    x = np.linspace(x_labels[0], x_labels[-1], len(x_labels) * 20)
    regression_curve1 = [coefficient1 * (_x ** exponent1) for _x in x]
    regression_curve2 = [coefficient2 * (_x ** exponent2) for _x in x]
    
    # Extract metadata (same for all files)
    algorithm1 = data_points1[0]["algorithm"]
    algorithm2 = data_points2[0]["algorithm"]
    point_count = data_points1[0]["point_count"]
    dimension_count1 = data_points1[0]["dimension"]
    dimension_count2 = data_points1[0]["dimension"]
    
    # Plot data
    plt.figure(figsize=(10, 5))
    plt.plot(x_labels, avg_times1, marker='d', label='Avg Time 1')
    plt.plot(x_labels, avg_times2, marker='d', label='Avg Time 2')

    plt.plot(x, regression_curve1, label=f"Curve y = {coefficient1:.2e} * x^{exponent1:.2f}")
    plt.plot(x, regression_curve2, label=f"Curve y = {coefficient2:.2e} * x^{exponent2:.2f}")
    
    # Labels and title

    plt.xlabel("Polyline Length")
    plt.ylabel("Time (seconds)")
    plt.title(f"Performance of {algorithm1} and {algorithm2}\nDimensions: {dimension_count1}, {dimension_count2}")
    plt.legend()
    ax = plt.gca()
    ax.grid(True)
    # ax.set_yscale('function', functions=(cubic_transform, inverse_cubic_transform))
    # time_ticks = np.arange(0, max(avg_times1[-1], avg_times2[-1]) + 0.21, 0.1)
    # time_ticks = [inverse_cubic_transform(t) for t in time_ticks]
    # ax.set_yticks(time_ticks)
    # ax.set_yscale('log')
    plt.grid()
    
    # Show plot
    plt.show()



if __name__ == "__main__":

    if len(sys.argv) == 2:
        plot_data(sys.argv[1])
    elif len(sys.argv) == 3:
        plot_data_compare(sys.argv[1], sys.argv[2])
    else:
        print("Usage: python measurement_plot.py <directory> or python measurement_plot.py <directory> <directory>")
        sys.exit(1)




