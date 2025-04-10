import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def read_csv_file(file_path):
    """Read a CSV file and return a DataFrame."""
    try:
        return pd.read_csv(file_path, header=None)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def compare_signals():
    """Compare the generated signals from Python and C versions."""
    c_signal = read_csv_file("signal.csv")
    py_signal = read_csv_file("python_signal.csv")
    
    if c_signal is None or py_signal is None:
        print("Signal files not found. Run both test programs first.")
        return
    
    # Check if signals have the same length
    if len(c_signal) != len(py_signal):
        print(f"Signal length mismatch: C={len(c_signal)}, Python={len(py_signal)}")
    
    # Calculate differences
    c_values = c_signal[1].values
    py_values = py_signal[1].values
    
    # Check for mismatches
    mismatches = np.where(c_values != py_values)[0]
    if len(mismatches) == 0:
        print("Signals match perfectly!")
    else:
        print(f"Found {len(mismatches)} signal value mismatches")
        for idx in mismatches[:10]:  # Show first 10 mismatches
            print(f"  Index {idx}: C={c_values[idx]}, Python={py_values[idx]}")
    
    # Plot signals for visual comparison
    plt.figure(figsize=(12, 6))
    plt.plot(c_values, label='C Signal')
    plt.plot(py_values, label='Python Signal', linestyle='--')
    plt.legend()
    plt.title('Signal Comparison')
    plt.savefig('signal_comparison.png')
    plt.close()

def compare_peak_results(test_name):
    """Compare peak detection results for a specific test."""
    c_file = f"peaks_{test_name}.csv"
    py_file = f"python_peaks_{test_name}.csv"
    
    if not os.path.exists(c_file) or not os.path.exists(py_file):
        print(f"Peak files for {test_name} test not found.")
        return
    
    c_results = read_csv_file(c_file)
    py_results = read_csv_file(py_file)
    
    if c_results is None or py_results is None:
        return
    
    # Get peak indices
    c_peaks = set(c_results[0].values)
    py_peaks = set(py_results[0].values)
    
    # Compare peak indices
    common_peaks = c_peaks.intersection(py_peaks)
    c_only = c_peaks - py_peaks
    py_only = py_peaks - c_peaks
    
    print(f"\nComparing {test_name} test results:")
    print(f"  Common peaks: {len(common_peaks)}")
    print(f"  C only peaks: {len(c_only)}")
    if len(c_only) > 0:
        print(f"    Indices: {sorted(c_only)}")
    print(f"  Python only peaks: {len(py_only)}")
    if len(py_only) > 0:
        print(f"    Indices: {sorted(py_only)}")
    
    # Calculate agreement percentage
    total_peaks = len(c_peaks.union(py_peaks))
    if total_peaks > 0:
        agreement = len(common_peaks) / total_peaks * 100
        print(f"  Agreement: {agreement:.2f}%")
    
    # Compare properties for common peaks
    if len(common_peaks) > 0 and c_results.shape[1] > 2 and py_results.shape[1] > 2:
        # Set up dictionaries for property comparison
        c_props = {c_results.iloc[i, 0]: c_results.iloc[i, 2:] for i in range(len(c_results))}
        py_props = {py_results.iloc[i, 0]: py_results.iloc[i, 2:] for i in range(len(py_results))}
        
        # Check property differences for common peaks
        prop_diffs = []
        for peak in common_peaks:
            c_peak_idx = c_results[0].values.tolist().index(peak)
            py_peak_idx = py_results[0].values.tolist().index(peak)
            
            c_props_row = c_results.iloc[c_peak_idx, 2:].values
            py_props_row = py_results.iloc[py_peak_idx, 2:].values
            
            # Match lengths if needed
            min_len = min(len(c_props_row), len(py_props_row))
            if min_len > 0:
                diff = np.abs(c_props_row[:min_len] - py_props_row[:min_len])
                prop_diffs.append((peak, diff.tolist()))
        
        # Report property differences
        if prop_diffs:
            print("  Property differences for common peaks:")
            for peak, diffs in prop_diffs[:5]:  # Show first 5
                print(f"    Peak at index {peak}, differences: {diffs}")

def run_comparison():
    """Run all comparisons."""
    print("Comparing C and Python find_peaks results")
    print("=========================================")
    
    # Compare signals
    print("\nComparing test signals:")
    compare_signals()
    
    # Compare peak results for each test
    for test in ["distance", "prominence", "width", "combined"]:
        compare_peak_results(test)

if __name__ == "__main__":
    run_comparison() 