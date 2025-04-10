#ifndef FIND_PEAKS_MTK_H
#define FIND_PEAKS_MTK_H

#include <stdbool.h>
#include <stdint.h>

// MTK specific type definitions
typedef int kal_int32;
typedef unsigned int kal_uint32;

/**
 * Peak results structure, containing found peaks and their properties
 */
typedef struct {
    int* indices;        // Indices of peaks
    int size;            // Number of peaks
    int capacity;        // Currently allocated capacity
    
    // Properties related to prominence
    uint32_t* prominences;
    int* left_bases;
    int* right_bases;
    
    // Properties related to width
    uint32_t* widths;
    uint32_t* width_heights;
    uint32_t* left_ips;
    uint32_t* right_ips;
    
    // Properties related to plateau size
    int* plateau_sizes;
    int* left_edges;
    int* right_edges;
} PeakResults;

/**
 * Find peaks in a signal
 *
 * @param x Input signal
 * @param size Size of the signal
 * @param distance Minimum distance between neighboring peaks, ignored if <= 1
 * @param min_prominence Minimum prominence of peaks, ignored if <= 0
 * @param max_prominence Maximum prominence of peaks, ignored if <= 0
 * @param min_width Minimum width of peaks, ignored if <= 0
 * @param max_width Maximum width of peaks, ignored if <= 0
 * @param min_plateau_size Minimum size of flat top of peaks, ignored if <= 0
 * @param max_plateau_size Maximum size of flat top of peaks, ignored if <= 0
 * @param wlen Window length for prominence calculation, uses full length if <= 0
 * @param rel_height Relative height for width calculation (0~1), default 0.5
 * @return Peak results structure
 */
PeakResults* find_peaks(const uint32_t* x, int size, 
                        int distance, 
                        uint32_t min_prominence, uint32_t max_prominence, 
                        uint32_t min_width, uint32_t max_width, 
                        int min_plateau_size, int max_plateau_size,
                        int wlen, uint32_t rel_height);

/**
 * Free PeakResults structure and all allocated resources
 */
void free_peak_results(PeakResults* results);

// MTK memory management functions declarations (for reference only)
extern void* get_ctrl_buffer(unsigned int size);
extern void free_ctrl_buffer(void* buff_ptr);

#endif /* FIND_PEAKS_MTK_H */ 