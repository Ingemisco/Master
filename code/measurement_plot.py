import os
import json
import matplotlib.pyplot as plt
from matplotlib.scale import FuncScale
import sys
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit

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

    
    # Extract metadata (same for all files)
    algorithm = data_points[0]["algorithm"]
    point_count = data_points[0]["point_count"]
    dimension_count = data_points[0]["dimension"]
    
    # Plot data
    plt.figure(figsize=(10, 5))
    plt.plot(x_labels, min_times, marker='o', label='Min Time')
    plt.plot(x_labels, max_times, marker='s', label='Max Time')
    plt.plot(x_labels, avg_times, marker='d', label='Avg Time')


    initial_guess = [1.0e-9, 3.5]
    params_opt, params_cov = curve_fit(
        lambda x, a, b: a * (x**b), 
        x_labels, 
        avg_times, 
        p0=initial_guess
    )
    coefficient, exponent = params_opt
    x = np.linspace(x_labels[0], x_labels[-1], 100)
    regression_curve = [coefficient * (_x ** exponent) for _x in x]
    plt.plot(x, regression_curve, label=f"Curve y = {coefficient:.2e} * x^{exponent:.2f}")
    
    # Labels and title

    plt.xlabel("Polyline Length")
    plt.ylabel("Time (seconds)")
    plt.title(f"Performance of {algorithm}\nDimensions: {dimension_count}")
    plt.legend()
    ax = plt.gca()
    ax.grid(True)
    plt.grid()
    
    # Show plot
    plt.show()

def joint_model(x, a1, a2, b):
    # Split input x into x1 and x2 (assuming equal lengths)
    n = len(x) // 2
    x1, x2 = x[:n], x[n:]
    
    # Compute predictions for both datasets
    y1_pred = a1 * (x1 ** b)
    y2_pred = a2 * (x2 ** b)
    
    # Concatenate predictions
    return np.concatenate([y1_pred, y2_pred])

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

    x_combined = np.concatenate([x_labels, x_labels])
    y_combined = np.concatenate([avg_times1, avg_times2])
    initial_guess = [1.0e-9, 1.0e-9, 3.5]
    params_opt, params_cov = curve_fit(
        joint_model, 
        x_combined, 
        y_combined, 
        p0=initial_guess
    )
    coefficient1, coefficient2, exponent = params_opt

    x = np.linspace(x_labels[0], x_labels[-1], len(x_labels) * 20)
    regression_curve1 = [coefficient1 * (_x ** exponent) for _x in x]
    regression_curve2 = [coefficient2 * (_x ** exponent) for _x in x]
    
    # Extract metadata (same for all files)
    algorithm1 = data_points1[0]["algorithm"]
    algorithm2 = data_points2[0]["algorithm"]
    point_count = data_points1[0]["point_count"]
    dimension_count1 = data_points1[0]["dimension"]
    dimension_count2 = data_points1[0]["dimension"]
    
    # Plot data
    plt.figure(figsize=(10, 5))
    plt.plot(x_labels, avg_times1, marker='d', label=f'Avg Time {directory1.split("/")[-1]}')
    plt.plot(x_labels, avg_times2, marker='d', label=f'Avg Time {directory2.split("/")[-1]}')

    plt.plot(x, regression_curve1, label=f"Curve y = {coefficient1:.2e} * x^{exponent:.2f}")
    plt.plot(x, regression_curve2, label=f"Curve y = {coefficient2:.2e} * x^{exponent:.2f}")
    
    # Labels and title

    plt.xlabel("Polyline Length")
    plt.ylabel("Time (seconds)")
    plt.title(f"Performance of {algorithm1} and {algorithm2}\nDimensions: {dimension_count1}, {dimension_count2}")
    plt.legend()
    ax = plt.gca()
    ax.grid(True)
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




