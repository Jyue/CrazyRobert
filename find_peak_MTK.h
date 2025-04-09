/**
 * find_peak.h - C implementation of peak detection functions
 * Converted from Python's scipy.signal module find_peaks functionality
 */

#ifndef FIND_PEAK_H
#define FIND_PEAK_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <string.h>

// Define the structure for returning peak properties
typedef struct {
    double* peak_heights;       // Peak heights
    double* prominences;        // Prominences
    size_t* left_bases;         // Left bases
    size_t* right_bases;        // Right bases
    double* widths;             // Widths
    double* width_heights;      // Width heights
    double* left_ips;           // Left intersection points
    double* right_ips;          // Right intersection points
    size_t* plateau_sizes;      // Plateau sizes
    size_t* left_edges;         // Left edges
    size_t* right_edges;        // Right edges
    size_t count;               // Number of peaks
} PeakProperties;

// Main peak finding function
/**
 * Find peaks in a signal
 * @param x Input signal array
 * @param n Signal length
 * @param height_min Minimum height (can be NULL)
 * @param distance Minimum distance between peaks (can be 0 to ignore)
 * @param prominence_min Minimum prominence (can be NULL)
 * @param width_min Minimum width (can be NULL)
 * @param wlen Window length for calculating prominence
 * @param rel_height Relative height for calculating width (between 0-1)
 * @param plateau_size_min Minimum plateau size (can be NULL)
 * @param peaks Output array of peak indices (must be freed using free)
 * @param properties Output peak properties (must be freed using free_peak_properties)
 * @return Number of peaks
 */
size_t find_peaks(const double* x, size_t n,
                 const double* height_min,
                 size_t distance,
                 const double* prominence_min,
                 const double* width_min,
                 size_t wlen, double rel_height,
                 const size_t* plateau_size_min,
                 size_t** peaks, PeakProperties* properties);

/**
 * Free memory in the peak properties structure
 * @param props The peak properties structure to free
 */
void free_peak_properties(PeakProperties* props);

/**
 * Calculate peak prominences
 * @param x Input signal array
 * @param n Signal length
 * @param peaks Peak indices array
 * @param n_peaks Number of peaks
 * @param wlen Window length (can be 0 to use the entire signal)
 * @param prominences Output prominences array (must be pre-allocated)
 * @param left_bases Output left base indices (must be pre-allocated)
 * @param right_bases Output right base indices (must be pre-allocated)
 */
void peak_prominences(const double* x, size_t n,
                     const size_t* peaks, size_t n_peaks,
                     size_t wlen,
                     double* prominences, 
                     size_t* left_bases, 
                     size_t* right_bases);

/**
 * Calculate peak widths
 * @param x Input signal array
 * @param n Signal length
 * @param peaks Peak indices array
 * @param n_peaks Number of peaks
 * @param rel_height Relative height (between 0-1)
 * @param prominences Peak prominences array
 * @param left_bases Left base indices array
 * @param right_bases Right base indices array
 * @param widths Output widths array (must be pre-allocated)
 * @param width_heights Output width heights array (must be pre-allocated)
 * @param left_ips Output left intersection points array (must be pre-allocated)
 * @param right_ips Output right intersection points array (must be pre-allocated)
 */
void peak_widths(const double* x, size_t n,
                const size_t* peaks, size_t n_peaks,
                double rel_height,
                const double* prominences,
                const size_t* left_bases,
                const size_t* right_bases,
                double* widths,
                double* width_heights,
                double* left_ips,
                double* right_ips);

#endif /* FIND_PEAK_H */ 