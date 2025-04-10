# C Implementation of Peak Finding Library

This library provides a C implementation of a function similar to Python SciPy's `find_peaks` function, used for finding peaks in one-dimensional signals.

## Features

Implements four main parameters of SciPy's `find_peaks` function:

1. `distance` - Minimum distance between neighboring peaks
2. `prominence` - Prominence constraints for peaks
3. `width` - Width constraints for peaks
4. `plateau_size` - Size constraints for flat tops of peaks

## Compiling and Running

Compile the test program using the provided Makefile:

```bash
make
```

Run the test program:

```bash
./test_peaks
```

After execution, several CSV files will be generated in the current directory, containing the test signal and peak information that can be viewed using spreadsheet software or plotting tools.

## Usage

1. Include the header file in your project:

```c
#include "find_peaks.h"
```

2. Call the `find_peaks` function:

```c
PeakResults* results = find_peaks(
    signal,     // Input signal
    signal_size,// Signal size
    distance,   // Minimum distance between peaks
    min_prominence, max_prominence, // Prominence range
    min_width, max_width,           // Width range
    min_plateau_size, max_plateau_size, // Plateau size range
    wlen,       // Window length for prominence calculation
    rel_height  // Relative height for width calculation
);
```

3. Use the results and free memory when done:

```c
// Use the results
for (int i = 0; i < results->size; i++) {
    int peak_index = results->indices[i];
    uint32_t peak_height = signal[peak_index];
    // Other processing...
}

// Free memory
free_peak_results(results);
```

## Parameter Description

- `distance`: Minimum distance between neighboring peaks, in samples. Ignored if <= 1.
- `min_prominence`, `max_prominence`: Minimum/maximum prominence of peaks. Ignored if <= 0.
- `min_width`, `max_width`: Minimum/maximum width of peaks, in samples. Ignored if <= 0.
- `min_plateau_size`, `max_plateau_size`: Minimum/maximum size of flat tops of peaks, in samples. Ignored if <= 0.
- `wlen`: Window length used for prominence calculation. Uses full length if <= 0.
- `rel_height`: Relative height used for width calculation (0~100), with 50 representing 0.5.

## Return Value

The function returns a `PeakResults` structure containing the found peaks and their properties:

- `indices`: Indices of the peaks
- `size`: Number of peaks found
- `prominences`: Prominences of the peaks
- `widths`: Widths of the peaks
- `plateau_sizes`: Sizes of flat tops of the peaks
- And other related properties

## Data Types

All signal values, prominences, widths and other floating-point values are represented using `uint32_t` type instead of floating point, providing better compatibility with embedded systems.

## Validation Against Python Version

To validate that this C implementation produces results consistent with the original SciPy Python implementation, follow these steps:

1. Ensure you have Python installed with NumPy, SciPy, Matplotlib, and Pandas packages.

2. Run the comparison script:

```bash
run_comparison.bat
```

This script will:
- Compile and run the C version of the code
- Run the Python reference implementation
- Compare the results between the two versions
- Generate comparison reports and visualizations

The comparison examines:
- Signal matching between C and Python versions
- Peak detection agreement for various parameter sets
- Consistency of calculated properties (prominence, width, etc.)

If there are significant differences, the script will highlight them to help identify implementation issues.

### Understanding Differences

Small differences between the C and Python implementations are expected due to:
- Different numeric precision (uint32_t vs floating-point)
- Implementation variations in interpolation and property calculations
- Platform-specific behavior

The most important validation metric is that both implementations identify the same peak locations with similar prominence and width values.
![image](https://github.com/user-attachments/assets/25b71f07-bedb-47bc-9898-a0fba10b4fb7)


## References

This implementation references SciPy's `find_peaks` function, but has been simplified to retain only the four most important filtering parameters. 
