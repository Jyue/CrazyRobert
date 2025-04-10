## MTK Adaptation

The peak finding library has been adapted for use with MTK's memory management system. This adaptation includes the following changes:

### Memory Management

1. All memory allocation and deallocation in the MTK version use the custom memory management functions `get_ctrl_buffer()` and `free_ctrl_buffer()` instead of standard C functions like `malloc()`, `calloc()`, and `free()`.

2. These functions serve as wrappers for memory management that can be integrated with MTK's internal memory management system.

### File Structure

The MTK version consists of three main files:

1. `find_peaks_MTK.h` - Header file containing type definitions, function declarations, and external declarations of MTK memory management functions.

2. `find_peaks_MTK.c` - Implementation file with all peak finding logic, using MTK memory management functions for all memory operations.

3. `test_peaks_MTK.c` - Test program that demonstrates using the peak finding functionality with the MTK memory management system.

### Test Implementation

The test implementation in `test_peaks_MTK.c` includes:

1. Temporary implementations of `get_ctrl_buffer()` and `free_ctrl_buffer()` that wrap standard C memory functions for testing purposes. These will be replaced by actual MTK functions in the production environment.

2. Functions to generate test signals, save signals to CSV files, and save peak results to CSV files.

3. A test suite that evaluates different peak finding parameters:
   - Distance between peaks
   - Peak prominence
   - Peak width
   - Plateau size
   - Combined parameters

### Usage in MTK Environment

To use this library in an MTK environment:

1. Include the header file in your project:
```c
#include "find_peaks_MTK.h"
```

2. Ensure that `get_ctrl_buffer()` and `free_ctrl_buffer()` are properly linked to your MTK system's memory management.

3. Call the peak finding function with appropriate parameters:
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

4. Process the results and free memory as needed:
```c
// Process results
for (int i = 0; i < results->size; i++) {
    // Access peak data
    int peak_index = results->indices[i];
    // ...
}

// Free memory
free_peak_results(results);
```
