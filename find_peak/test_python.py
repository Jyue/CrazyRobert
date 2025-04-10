import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import csv

# Generate the same test signal as in C version
def generate_test_signal(size=100):
    x = np.arange(size)
    signal = np.zeros(size, dtype=np.uint32)
    
    # Add Gaussian peaks (same as C version)
    signal += np.uint32(5 * np.exp(-0.1 * (x - 20) ** 2))
    signal += np.uint32(3 * np.exp(-0.1 * (x - 50) ** 2))
    signal += np.uint32(7 * np.exp(-0.1 * (x - 80) ** 2))
    
    # Add plateau peak (same as C version)
    signal[30:36] = 4
    
    # Add noise (approximating the C version noise)
    np.random.seed(42)  # For reproducibility
    noise = np.uint32((np.random.random(size) - 0.5) * 0.5)
    signal += noise
    
    return signal

# Save signal to CSV (same format as C version)
def save_signal_to_csv(signal, filename):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for i, val in enumerate(signal):
            writer.writerow([i, val])

# Save peaks to CSV (same format as C version)
def save_peaks_to_csv(signal, peak_indices, properties, filename):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for i, idx in enumerate(peak_indices):
            row = [idx, signal[idx]]
            
            # Add properties if available
            if 'prominences' in properties:
                row.append(int(properties['prominences'][i]))
            if 'widths' in properties:
                row.append(int(properties['widths'][i]))
            if 'plateau_sizes' in properties:
                row.append(int(properties.get('plateau_sizes', [0]*len(peak_indices))[i]))
                
            writer.writerow(row)

# Find plateaus (flat tops) in signal
def find_plateaus(signal, min_plateau_size=1):
    plateaus = []
    plateau_sizes = []
    
    i = 0
    while i < len(signal) - 1:
        # Check for start of plateau
        if i > 0 and signal[i] > signal[i-1]:
            # Found potential start
            start_idx = i
            current_val = signal[i]
            
            # Find plateau end
            j = i + 1
            while j < len(signal) and signal[j] == current_val:
                j += 1
            
            plateau_size = j - start_idx
            
            # Check if end of plateau is a peak (value decreases after)
            if j < len(signal) and signal[j] < current_val and plateau_size >= min_plateau_size:
                # Midpoint of plateau
                mid_idx = start_idx + plateau_size // 2
                plateaus.append(mid_idx)
                plateau_sizes.append(plateau_size)
                
            i = j  # Skip to end of plateau
        else:
            i += 1
    
    return np.array(plateaus), {'plateau_sizes': np.array(plateau_sizes)}

# Run tests with the same parameters as C version
def run_tests():
    signal = generate_test_signal(100)
    
    # Save the signal to CSV for comparison
    save_signal_to_csv(signal, "python_signal.csv")
    
    # Test distance parameter
    peaks1, properties1 = find_peaks(signal, distance=15)
    print(f"Testing distance parameter:")
    print(f"Found {len(peaks1)} peaks using distance=15")
    for i, idx in enumerate(peaks1):
        print(f"  Peak #{i+1}: index={idx}, height={signal[idx]}")
    save_peaks_to_csv(signal, peaks1, properties1, "python_peaks_distance.csv")
    
    # Test prominence parameter
    peaks2, properties2 = find_peaks(signal, prominence=3)
    print(f"\nTesting prominence parameter:")
    print(f"Found {len(peaks2)} peaks using min_prominence=3")
    for i, idx in enumerate(peaks2):
        print(f"  Peak #{i+1}: index={idx}, height={signal[idx]}, prominence={int(properties2['prominences'][i])}")
    save_peaks_to_csv(signal, peaks2, properties2, "python_peaks_prominence.csv")
    
    # Test width parameter
    peaks3, properties3 = find_peaks(signal, width=5)
    print(f"\nTesting width parameter:")
    print(f"Found {len(peaks3)} peaks using min_width=5")
    for i, idx in enumerate(peaks3):
        print(f"  Peak #{i+1}: index={idx}, height={signal[idx]}, width={int(properties3['widths'][i])}")
    save_peaks_to_csv(signal, peaks3, properties3, "python_peaks_width.csv")
    
    # Test plateau_size parameter - SciPy doesn't have this, so we'll implement it
    peaks4, properties4 = find_plateaus(signal, min_plateau_size=3)
    print(f"\nTesting plateau size parameter:")
    print(f"Found {len(peaks4)} peaks using min_plateau_size=3")
    for i, idx in enumerate(peaks4):
        print(f"  Peak #{i+1}: index={idx}, height={signal[idx]}, plateau_size={properties4['plateau_sizes'][i]}")
    save_peaks_to_csv(signal, peaks4, properties4, "python_peaks_plateau.csv")
    
    # Test combined parameters
    peaks5, properties5 = find_peaks(signal, distance=10, prominence=2, width=3)
    print(f"\nTesting combined parameters:")
    print(f"Found {len(peaks5)} peaks using combined parameters")
    for i, idx in enumerate(peaks5):
        print(f"  Peak #{i+1}: index={idx}, height={signal[idx]}, prominence={int(properties5['prominences'][i])}, width={int(properties5['widths'][i])}")
    save_peaks_to_csv(signal, peaks5, properties5, "python_peaks_combined.csv")
    
    # Create a plot to visualize
    plt.figure(figsize=(10, 6))
    plt.plot(signal, label='Signal')
    plt.plot(peaks5, signal[peaks5], 'x', color='red', label='Detected Peaks')
    plt.legend()
    plt.title('Python find_peaks result')
    plt.savefig('python_peaks_plot.png')
    plt.close()

if __name__ == "__main__":
    run_tests() 