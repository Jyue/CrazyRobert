import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths

# Load test signal
data = np.loadtxt('test_signal.csv', delimiter=',', skiprows=1, usecols=(1,))

# Use scipy.signal.find_peaks to find peaks
peaks, properties = find_peaks(data, height=0, distance=50, prominence=0.1)

# Calculate peak widths
widths, width_heights, left_ips, right_ips = peak_widths(data, peaks, rel_height=0.5)

# Save results
with open('python_find_peaks_results.csv', 'w') as f:
    f.write('index,signal,is_peak,prominence,width\n')
    for i in range(len(data)):
        is_peak = 1 if i in peaks else 0
        row = f'{i},{data[i]},{is_peak}'
        if is_peak:
            peak_idx = np.where(peaks == i)[0][0]
            row += f',{properties["prominences"][peak_idx]},{widths[peak_idx]}'
        else:
            row += ',,'
        f.write(row + '\n')

# Compare results
plt.figure(figsize=(12, 6))
plt.plot(data)
plt.plot(peaks, data[peaks], 'ro', label='Python Peaks')

# Load C version results
try:
    c_results = np.loadtxt('c_find_peaks_results.csv', delimiter=',', skiprows=1, usecols=(0,1,2))
    c_peaks = c_results[c_results[:, 2] == 1, 0].astype(int)
    plt.plot(c_peaks, data[c_peaks], 'gx', label='C Peaks')
except Exception as e:
    print(f"Could not load C results: {e}")
    c_peaks = []

plt.legend()
plt.title('Python vs C find_peaks Comparison')
plt.savefig('peaks_comparison.png')
plt.show()

# Print comparison results
print(f'Python found {len(peaks)} peaks')
print(f'C found {len(c_peaks)} peaks')
print(f'Number of common peaks: {len(np.intersect1d(peaks, c_peaks))}')
