import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read signal data
def read_signal_data(filename):
    try:
        df = pd.read_csv(filename, header=None)
        signal = df[1].values
        return signal
    except Exception as e:
        print(f"Error reading signal file {filename}: {e}")
        return None

# Read peak data
def read_peak_data(filename):
    try:
        df = pd.read_csv(filename, header=None)
        peak_indices = df[0].values
        return peak_indices
    except Exception as e:
        print(f"Error reading peak file {filename}: {e}")
        return None

# Plot peak comparison on the same chart
def plot_comparison(signal, c_peaks, py_peaks, test_name, output_filename):
    plt.figure(figsize=(12, 6))
    
    # Plot signal
    plt.plot(signal, label='Signal')
    
    # Plot C version peaks
    if c_peaks is not None:
        plt.plot(c_peaks, signal[c_peaks], 'o', markersize=10, 
                 markerfacecolor='none', markeredgecolor='r', markeredgewidth=2,
                 label='C Peaks')
    
    # Plot Python version peaks
    if py_peaks is not None:
        plt.plot(py_peaks, signal[py_peaks], 'x', markersize=8, 
                 color='b', markeredgewidth=2, 
                 label='Python Peaks')
    
    plt.title(f'{test_name} Parameter Test - Peak Comparison')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlabel('Signal Index')
    plt.ylabel('Signal Value')
    
    # Label each peak index
    all_peaks = set()
    if c_peaks is not None:
        all_peaks.update(c_peaks)
    if py_peaks is not None:
        all_peaks.update(py_peaks)
    
    for peak in all_peaks:
        plt.annotate(f'{peak}', 
                    xy=(peak, signal[peak]), 
                    xytext=(0, 10),
                    textcoords='offset points',
                    ha='center',
                    fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_filename)
    print(f"Generated comparison chart: {output_filename}")

def main():
    # Read signal
    signal = read_signal_data('signal.csv')
    if signal is None:
        return
    
    # Test names and corresponding files
    tests = [
        ("Distance", "peaks_distance.csv", "python_peaks_distance.csv", "distance_comparison.png"),
        ("Prominence", "peaks_prominence.csv", "python_peaks_prominence.csv", "prominence_comparison.png"),
        ("Width", "peaks_width.csv", "python_peaks_width.csv", "width_comparison.png"),
        ("Plateau", "peaks_plateau.csv", "python_peaks_plateau.csv", "plateau_comparison.png"),
        ("Combined", "peaks_combined.csv", "python_peaks_combined.csv", "combined_comparison.png")
    ]
    
    # Generate comparison charts for each test
    for test_name, c_file, py_file, output_file in tests:
        c_peaks = read_peak_data(c_file)
        py_peaks = read_peak_data(py_file)
        
        if c_peaks is not None or py_peaks is not None:
            plot_comparison(signal, c_peaks, py_peaks, test_name, output_file)
    
    # Generate a summary chart with all test results
    plt.figure(figsize=(15, 10))
    
    # Create subplots
    for i, (test_name, c_file, py_file, _) in enumerate(tests, 1):
        plt.subplot(3, 2, i)
        
        c_peaks = read_peak_data(c_file)
        py_peaks = read_peak_data(py_file)
        
        # Plot signal
        plt.plot(signal, label='Signal')
        
        # Plot C version peaks
        if c_peaks is not None:
            plt.plot(c_peaks, signal[c_peaks], 'o', markersize=8, 
                     markerfacecolor='none', markeredgecolor='r', markeredgewidth=2,
                     label='C Peaks')
        
        # Plot Python version peaks
        if py_peaks is not None:
            plt.plot(py_peaks, signal[py_peaks], 'x', markersize=6, 
                     color='b', markeredgewidth=2, 
                     label='Python Peaks')
        
        plt.title(f'{test_name} Parameter Test')
        plt.legend(loc='upper right', fontsize=8)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        if i > 3:  # Add x-axis label only to bottom subplots
            plt.xlabel('Signal Index')
        
        if i % 2 == 1:  # Add y-axis label only to left subplots
            plt.ylabel('Signal Value')
    
    plt.tight_layout()
    plt.savefig('all_tests_comparison.png')
    print("Generated summary chart: all_tests_comparison.png")

if __name__ == "__main__":
    main() 