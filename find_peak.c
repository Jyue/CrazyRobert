/**
 * find_peak.c - Implementation of C peak detection functions
 * Converted from Python's scipy.signal module find_peaks functionality
 */

#include "find_peak.h"

// Initialize peak properties structure
static void init_peak_properties(PeakProperties* props, size_t count) {
    props->count = count;
    props->peak_heights = NULL;
    props->prominences = NULL;
    props->left_bases = NULL;
    props->right_bases = NULL;
    props->widths = NULL;
    props->width_heights = NULL;
    props->left_ips = NULL;
    props->right_ips = NULL;
    props->plateau_sizes = NULL;
    props->left_edges = NULL;
    props->right_edges = NULL;
}

// Free memory in the peak properties structure
void free_peak_properties(PeakProperties* props) {
    if (props == NULL) return;
    
    free(props->peak_heights);
    free(props->prominences);
    free(props->left_bases);
    free(props->right_bases);
    free(props->widths);
    free(props->width_heights);
    free(props->left_ips);
    free(props->right_ips);
    free(props->plateau_sizes);
    free(props->left_edges);
    free(props->right_edges);
    
    // Reset all pointers
    init_peak_properties(props, 0);
}

// Sort array and return sorted indices
static void argsort(const double* array, size_t n, size_t* indices) {
    // Initialize index array
    for (size_t i = 0; i < n; i++) {
        indices[i] = i;
    }
    
    // Use insertion sort for indices
    for (size_t i = 1; i < n; i++) {
        size_t j = i;
        while (j > 0 && array[indices[j-1]] > array[indices[j]]) {
            // Swap indices
            size_t temp = indices[j];
            indices[j] = indices[j-1];
            indices[j-1] = temp;
            j--;
        }
    }
}

// Find local maxima
static size_t _local_maxima_1d(const double* x, size_t n, 
                              size_t** peaks, 
                              size_t** left_edges, 
                              size_t** right_edges) {
    // Allocate temporary array to mark local maxima
    bool* is_max = (bool*)calloc(n, sizeof(bool));
    if (!is_max) {
        return 0; // Memory allocation failed
    }
    
    // Find all local maxima
    size_t max_count = 0;
    
    // Special handling for the first point
    if (n > 1 && x[0] > x[1]) {
        is_max[0] = true;
        max_count++;
    }
    
    // Process middle points
    for (size_t i = 1; i < n - 1; i++) {
        if (x[i] > x[i-1] && x[i] >= x[i+1]) {
            is_max[i] = true;
            max_count++;
        }
    }
    
    // Special handling for the last point
    if (n > 1 && x[n-1] > x[n-2]) {
        is_max[n-1] = true;
        max_count++;
    }
    
    // Allocate result arrays
    *peaks = (size_t*)malloc(max_count * sizeof(size_t));
    *left_edges = (size_t*)malloc(max_count * sizeof(size_t));
    *right_edges = (size_t*)malloc(max_count * sizeof(size_t));
    
    if (!*peaks || !*left_edges || !*right_edges) {
        free(is_max);
        free(*peaks);
        free(*left_edges);
        free(*right_edges);
        *peaks = NULL;
        *left_edges = NULL;
        *right_edges = NULL;
        return 0; // Memory allocation failed
    }
    
    // Fill result arrays
    size_t count = 0;
    
    // Iterate to find plateau edges
    for (size_t i = 0; i < n; i++) {
        if (is_max[i]) {
            // Find left edge
            size_t left = i;
            while (left > 0 && x[left-1] == x[i]) {
                left--;
            }
            
            // Find right edge
            size_t right = i;
            while (right < n - 1 && x[right+1] == x[i]) {
                right++;
            }
            
            // Calculate peak position (middle of plateau for plateaus)
            size_t peak = left + (right - left) / 2;
            
            (*peaks)[count] = peak;
            (*left_edges)[count] = left;
            (*right_edges)[count] = right;
            count++;
            
            // Skip over processed plateau
            i = right;
        }
    }
    
    free(is_max);
    return count;
}

// Calculate percentile of data
static double scoreatpercentile(const double* data, size_t n, double per) {
    if (n == 0) return 0.0;
    if (n == 1) return data[0];
    
    // Copy data for sorting
    double* sorted = (double*)malloc(n * sizeof(double));
    if (!sorted) return 0.0;
    
    memcpy(sorted, data, n * sizeof(double));
    
    // Simple insertion sort
    for (size_t i = 1; i < n; i++) {
        double key = sorted[i];
        int j = i - 1;
        while (j >= 0 && sorted[j] > key) {
            sorted[j + 1] = sorted[j];
            j--;
        }
        sorted[j + 1] = key;
    }
    
    // Calculate index
    double idx = per / 100.0 * (n - 1);
    size_t idx_floor = (size_t)floor(idx);
    size_t idx_ceil = (size_t)ceil(idx);
    
    double result;
    if (idx_floor == idx_ceil) {
        result = sorted[idx_floor];
    } else {
        double frac = idx - idx_floor;
        result = (1 - frac) * sorted[idx_floor] + frac * sorted[idx_ceil];
    }
    
    free(sorted);
    return result;
}

// Select peaks based on property (can handle double and size_t types)
static bool* _select_by_property_double(const double* property, size_t n,
                                    const double* pmin, const double* pmax) {
    bool* keep = (bool*)malloc(n * sizeof(bool));
    if (!keep) return NULL;
    
    // Initialize all to keep
    for (size_t i = 0; i < n; i++) {
        keep[i] = true;
    }
    
    // Apply minimum threshold
    if (pmin != NULL) {
        for (size_t i = 0; i < n; i++) {
            if (property[i] < *pmin) {
                keep[i] = false;
            }
        }
    }
    
    // Apply maximum threshold
    if (pmax != NULL) {
        for (size_t i = 0; i < n; i++) {
            if (property[i] > *pmax) {
                keep[i] = false;
            }
        }
    }
    
    return keep;
}

// Use for processing size_t type property
static bool* _select_by_property_size_t(const size_t* property, size_t n,
                                     const size_t* pmin, const size_t* pmax) {
    bool* keep = (bool*)malloc(n * sizeof(bool));
    if (!keep) return NULL;
    
    // Initialize all to keep
    for (size_t i = 0; i < n; i++) {
        keep[i] = true;
    }
    
    // Apply minimum threshold
    if (pmin != NULL) {
        for (size_t i = 0; i < n; i++) {
            if (property[i] < *pmin) {
                keep[i] = false;
            }
        }
    }
    
    // Apply maximum threshold
    if (pmax != NULL) {
        for (size_t i = 0; i < n; i++) {
            if (property[i] > *pmax) {
                keep[i] = false;
            }
        }
    }
    
    return keep;
}

// Select peaks based on peak threshold
static bool* _select_by_peak_threshold(const double* x, const size_t* peaks, size_t n,
                                     const double* tmin, const double* tmax,
                                     double** left_thresholds, double** right_thresholds) {
    *left_thresholds = (double*)malloc(n * sizeof(double));
    *right_thresholds = (double*)malloc(n * sizeof(double));
    bool* keep = (bool*)malloc(n * sizeof(bool));
    
    if (!*left_thresholds || !*right_thresholds || !keep) {
        free(*left_thresholds);
        free(*right_thresholds);
        free(keep);
        *left_thresholds = NULL;
        *right_thresholds = NULL;
        return NULL;
    }
    
    // Initialize all to keep
    for (size_t i = 0; i < n; i++) {
        keep[i] = true;
    }
    
    // Calculate left and right thresholds for each peak
    for (size_t i = 0; i < n; i++) {
        size_t idx = peaks[i];
        double left_threshold = 0, right_threshold = 0;
        
        // Ensure index is within valid range
        if (idx > 0) {
            left_threshold = x[idx] - x[idx-1];
        }
        
        if (idx < n - 1) {
            right_threshold = x[idx] - x[idx+1];
        }
        
        (*left_thresholds)[i] = left_threshold;
        (*right_thresholds)[i] = right_threshold;
        
        // Calculate minimum threshold
        double min_threshold = (left_threshold < right_threshold) ? left_threshold : right_threshold;
        
        // Calculate maximum threshold
        double max_threshold = (left_threshold > right_threshold) ? left_threshold : right_threshold;
        
        // Apply minimum threshold condition
        if (tmin != NULL && min_threshold < *tmin) {
            keep[i] = false;
        }
        
        // Apply maximum threshold condition
        if (tmax != NULL && max_threshold > *tmax) {
            keep[i] = false;
        }
    }
    
    return keep;
}

// Select peaks based on peak distance
static bool* _select_by_peak_distance(const size_t* peaks, size_t n, const double* values, size_t distance) {
    if (distance <= 1 || n == 0) {
        // If no distance limit or no peaks, keep all peaks
        bool* keep = (bool*)malloc(n * sizeof(bool));
        if (!keep) return NULL;
        
        for (size_t i = 0; i < n; i++) {
            keep[i] = true;
        }
        return keep;
    }
    
    // First sort peaks by peak height
    size_t* sorted_indices = (size_t*)malloc(n * sizeof(size_t));
    if (!sorted_indices) return NULL;
    
    argsort(values, n, sorted_indices);
    
    // Create result array, initialize all to keep
    bool* keep = (bool*)malloc(n * sizeof(bool));
    if (!keep) {
        free(sorted_indices);
        return NULL;
    }
    
    for (size_t i = 0; i < n; i++) {
        keep[i] = true;
    }
    
    // Iterate from high to low by peak height
    for (size_t i = n; i > 0; i--) {
        size_t idx = sorted_indices[i-1];
        
        if (keep[idx]) {
            // Mark peaks too close to be kept
            size_t peak_idx = peaks[idx];
            
            for (size_t j = 0; j < n; j++) {
                if (j != idx && keep[j]) {
                    size_t other_idx = peaks[j];
                    if (other_idx >= peak_idx && other_idx - peak_idx < distance) {
                        keep[j] = false;
                    } else if (peak_idx >= other_idx && peak_idx - other_idx < distance) {
                        keep[j] = false;
                    }
                }
            }
        }
    }
    
    free(sorted_indices);
    return keep;
}

// Calculate peak prominence
void peak_prominences(const double* x, size_t n,
                     const size_t* peaks, size_t n_peaks,
                     size_t wlen,
                     double* prominences, 
                     size_t* left_bases, 
                     size_t* right_bases) {
    // If no peaks, return immediately
    if (n_peaks == 0) return;
    
    // If window size not specified, use entire signal
    size_t effective_wlen;
    if (wlen == 0 || wlen > n) {
        effective_wlen = n;
    } else {
        // Ensure window size is odd
        effective_wlen = wlen;
        if (effective_wlen % 2 == 0) {
            effective_wlen += 1;
        }
    }
    
    // Process each peak
    for (size_t i = 0; i < n_peaks; i++) {
        size_t peak_idx = peaks[i];
        double peak_height = x[peak_idx];
        
        // Define window range
        size_t window_start, window_end;
        
        if (effective_wlen < n) {
            // Use specified window size
            size_t half_wlen = effective_wlen / 2;
            window_start = (peak_idx > half_wlen) ? peak_idx - half_wlen : 0;
            window_end = (peak_idx + half_wlen < n) ? peak_idx + half_wlen : n - 1;
        } else {
            // Use entire signal
            window_start = 0;
            window_end = n - 1;
        }
        
        // Find minimum value in left side
        size_t left_min_idx = peak_idx;
        double left_min = peak_height;
        bool found_left_crossing = false;
        
        for (size_t j = peak_idx; j > window_start; j--) {
            // If higher peak found, stop searching
            if (x[j] > peak_height) {
                found_left_crossing = true;
                break;
            }
            
            // Update minimum value
            if (x[j] < left_min) {
                left_min = x[j];
                left_min_idx = j;
            }
        }
        
        // If window too small, check window start
        if (!found_left_crossing && window_start > 0) {
            if (x[window_start] < left_min) {
                left_min = x[window_start];
                left_min_idx = window_start;
            }
        }
        
        // Find minimum value in right side
        size_t right_min_idx = peak_idx;
        double right_min = peak_height;
        bool found_right_crossing = false;
        
        for (size_t j = peak_idx; j < window_end; j++) {
            // If higher peak found, stop searching
            if (x[j] > peak_height) {
                found_right_crossing = true;
                break;
            }
            
            // Update minimum value
            if (x[j] < right_min) {
                right_min = x[j];
                right_min_idx = j;
            }
        }
        
        // If window too small, check window end
        if (!found_right_crossing && window_end < n - 1) {
            if (x[window_end] < right_min) {
                right_min = x[window_end];
                right_min_idx = window_end;
            }
        }
        
        // Determine peak base height and prominence
        double base_height = (left_min > right_min) ? left_min : right_min;
        prominences[i] = peak_height - base_height;
        
        // Save base position
        left_bases[i] = left_min_idx;
        right_bases[i] = right_min_idx;
    }
}

// Linear interpolation to find intersection
static double _interpolate_intersection(double x1, double y1, double x2, double y2, double target_y) {
    // If two points have same y value, cannot interpolate, return one point
    if (y1 == y2) {
        return x1;
    }
    
    // Linear interpolation to calculate intersection
    return x1 + (x2 - x1) * (target_y - y1) / (y2 - y1);
}

// Calculate peak width
void peak_widths(const double* x, size_t n,
                const size_t* peaks, size_t n_peaks,
                double rel_height,
                const double* prominences,
                const size_t* left_bases __attribute__((unused)),
                const size_t* right_bases __attribute__((unused)),
                double* widths,
                double* width_heights,
                double* left_ips,
                double* right_ips) {
    // If no peaks, return immediately
    if (n_peaks == 0) return;
    
    // Calculate width for each peak
    for (size_t i = 0; i < n_peaks; i++) {
        size_t peak_idx = peaks[i];
        double prominence = prominences[i];
        double peak_height = x[peak_idx];
        
        // Calculate evaluation height
        double evaluation_height = peak_height - prominence * rel_height;
        width_heights[i] = evaluation_height;
        
        // Find intersection point
        size_t j_left = peak_idx;
        size_t j_right = peak_idx;
        
        // Search left for intersection
        while (j_left > 0 && x[j_left] > evaluation_height) {
            j_left--;
        }
        
        // Search right for intersection
        while (j_right < n - 1 && x[j_right] > evaluation_height) {
            j_right++;
        }
        
        // If reached signal boundary, use boundary value
        if (j_left == 0 && x[j_left] > evaluation_height) {
            left_ips[i] = 0.0;
        } else {
            // Calculate precise intersection point through linear interpolation
            size_t idx_left = j_left;
            size_t idx_right = j_left + 1;
            
            // Ensure not out of bounds
            if (idx_right >= n) {
                idx_right = n - 1;
            }
            
            left_ips[i] = _interpolate_intersection(idx_left, x[idx_left], idx_right, x[idx_right], evaluation_height);
        }
        
        // If reached signal boundary, use boundary value
        if (j_right == n - 1 && x[j_right] > evaluation_height) {
            right_ips[i] = n - 1.0;
        } else {
            // Calculate precise intersection point through linear interpolation
            size_t idx_left = j_right - 1;
            size_t idx_right = j_right;
            
            // Ensure not out of bounds
            if (idx_left >= n) {
                idx_left = n - 1;
            }
            
            right_ips[i] = _interpolate_intersection(idx_left, x[idx_left], idx_right, x[idx_right], evaluation_height);
        }
        
        // Calculate width
        widths[i] = right_ips[i] - left_ips[i];
    }
}

// Main peak finding function
size_t find_peaks(const double* x, size_t n,
                 const double* height_min,
                 size_t distance,
                 const double* prominence_min,
                 const double* width_min,
                 size_t wlen, double rel_height,
                 const size_t* plateau_size_min,
                 size_t** peaks, PeakProperties* properties) {
    // Initialize output
    *peaks = NULL;
    init_peak_properties(properties, 0);
    
    if (n == 0) {
        return 0;
    }
    
    // Check distance parameter
    if (distance == 0) {
        // Distance of 0 means no distance limit, which is allowed
    } else if (distance < 1) {
        fprintf(stderr, "Error: distance must be greater or equal to 1\n");
        return 0;
    }
    
    // Find all local maxima
    size_t* all_peaks = NULL;
    size_t* left_edges = NULL;
    size_t* right_edges = NULL;
    size_t peak_count = _local_maxima_1d(x, n, &all_peaks, &left_edges, &right_edges);
    
    if (peak_count == 0) {
        return 0; // No peaks found
    }
    
    // Save copy of all peaks for subsequent filtering
    size_t* current_peaks = (size_t*)malloc(peak_count * sizeof(size_t));
    if (!current_peaks) {
        free(all_peaks);
        free(left_edges);
        free(right_edges);
        return 0;
    }
    
    memcpy(current_peaks, all_peaks, peak_count * sizeof(size_t));
    size_t current_count = peak_count;
    
    // Evaluate plateau size
    if (plateau_size_min != NULL) {
        // Calculate plateau size
        size_t* plateau_sizes = (size_t*)malloc(current_count * sizeof(size_t));
        if (!plateau_sizes) {
            free(all_peaks);
            free(left_edges);
            free(right_edges);
            free(current_peaks);
            return 0;
        }
        
        for (size_t i = 0; i < current_count; i++) {
            plateau_sizes[i] = right_edges[i] - left_edges[i] + 1;
        }
        
        // Filter peaks that meet plateau size condition
        bool* keep = _select_by_property_size_t(plateau_sizes, current_count, 
                                                plateau_size_min, NULL);
        if (!keep) {
            free(all_peaks);
            free(left_edges);
            free(right_edges);
            free(current_peaks);
            free(plateau_sizes);
            return 0;
        }
        
        // Update peak list
        size_t new_count = 0;
        for (size_t i = 0; i < current_count; i++) {
            if (keep[i]) {
                current_peaks[new_count] = current_peaks[i];
                left_edges[new_count] = left_edges[i];
                right_edges[new_count] = right_edges[i];
                plateau_sizes[new_count] = plateau_sizes[i];
                new_count++;
            }
        }
        
        current_count = new_count;
        
        // Save plateau-related properties
        if (current_count > 0) {
            properties->plateau_sizes = (size_t*)malloc(current_count * sizeof(size_t));
            properties->left_edges = (size_t*)malloc(current_count * sizeof(size_t));
            properties->right_edges = (size_t*)malloc(current_count * sizeof(size_t));
            
            if (!properties->plateau_sizes || !properties->left_edges || !properties->right_edges) {
                free(all_peaks);
                free(left_edges);
                free(right_edges);
                free(current_peaks);
                free(plateau_sizes);
                free(keep);
                free_peak_properties(properties);
                return 0;
            }
            
            memcpy(properties->plateau_sizes, plateau_sizes, current_count * sizeof(size_t));
            memcpy(properties->left_edges, left_edges, current_count * sizeof(size_t));
            memcpy(properties->right_edges, right_edges, current_count * sizeof(size_t));
        }
        
        free(keep);
        free(plateau_sizes);
    }
    
    // Evaluate height condition
    if (height_min != NULL) {
        // Get peak heights
        double* peak_heights = (double*)malloc(current_count * sizeof(double));
        if (!peak_heights) {
            free(all_peaks);
            free(left_edges);
            free(right_edges);
            free(current_peaks);
            free_peak_properties(properties);
            return 0;
        }
        
        for (size_t i = 0; i < current_count; i++) {
            peak_heights[i] = x[current_peaks[i]];
        }
        
        // Filter peaks that meet height condition
        bool* keep = _select_by_property_double(peak_heights, current_count, height_min, NULL);
        if (!keep) {
            free(all_peaks);
            free(left_edges);
            free(right_edges);
            free(current_peaks);
            free(peak_heights);
            free_peak_properties(properties);
            return 0;
        }
        
        // Update peak list
        size_t new_count = 0;
        for (size_t i = 0; i < current_count; i++) {
            if (keep[i]) {
                current_peaks[new_count] = current_peaks[i];
                peak_heights[new_count] = peak_heights[i];
                if (properties->plateau_sizes) {
                    properties->plateau_sizes[new_count] = properties->plateau_sizes[i];
                    properties->left_edges[new_count] = properties->left_edges[i];
                    properties->right_edges[new_count] = properties->right_edges[i];
                }
                new_count++;
            }
        }
        
        current_count = new_count;
        
        // Save height properties
        if (current_count > 0) {
            properties->peak_heights = (double*)malloc(current_count * sizeof(double));
            if (!properties->peak_heights) {
                free(all_peaks);
                free(left_edges);
                free(right_edges);
                free(current_peaks);
                free(peak_heights);
                free(keep);
                free_peak_properties(properties);
                return 0;
            }
            
            memcpy(properties->peak_heights, peak_heights, current_count * sizeof(double));
        }
        
        free(keep);
        free(peak_heights);
    }
    
    // Evaluate distance condition
    if (distance > 1) {
        // Get peak heights (if not already calculated)
        double* peak_heights = NULL;
        if (properties->peak_heights) {
            peak_heights = properties->peak_heights;
        } else {
            peak_heights = (double*)malloc(current_count * sizeof(double));
            if (!peak_heights) {
                free(all_peaks);
                free(left_edges);
                free(right_edges);
                free(current_peaks);
                free_peak_properties(properties);
                return 0;
            }
            
            for (size_t i = 0; i < current_count; i++) {
                peak_heights[i] = x[current_peaks[i]];
            }
        }
        
        // Filter peaks that meet distance condition
        bool* keep = _select_by_peak_distance(current_peaks, current_count, peak_heights, distance);
        if (!keep) {
            if (!properties->peak_heights) {
                free(peak_heights);
            }
            free(all_peaks);
            free(left_edges);
            free(right_edges);
            free(current_peaks);
            free_peak_properties(properties);
            return 0;
        }
        
        // Update peak list
        size_t new_count = 0;
        for (size_t i = 0; i < current_count; i++) {
            if (keep[i]) {
                current_peaks[new_count] = current_peaks[i];
                if (properties->peak_heights) {
                    properties->peak_heights[new_count] = properties->peak_heights[i];
                }
                if (properties->plateau_sizes) {
                    properties->plateau_sizes[new_count] = properties->plateau_sizes[i];
                    properties->left_edges[new_count] = properties->left_edges[i];
                    properties->right_edges[new_count] = properties->right_edges[i];
                }
                new_count++;
            }
        }
        
        current_count = new_count;
        
        // If new peak_heights were created, free them
        if (!properties->peak_heights) {
            free(peak_heights);
        }
        
        free(keep);
    }
    
    // Evaluate prominence condition
    if (prominence_min != NULL || width_min != NULL) {
        // Calculate prominence
        double* prominences = (double*)malloc(current_count * sizeof(double));
        size_t* left_bases = (size_t*)malloc(current_count * sizeof(size_t));
        size_t* right_bases = (size_t*)malloc(current_count * sizeof(size_t));
        
        if (!prominences || !left_bases || !right_bases) {
            free(all_peaks);
            free(left_edges);
            free(right_edges);
            free(current_peaks);
            free(prominences);
            free(left_bases);
            free(right_bases);
            free_peak_properties(properties);
            return 0;
        }
        
        peak_prominences(x, n, current_peaks, current_count, wlen, prominences, left_bases, right_bases);
        
        // Filter peaks that meet prominence condition
        if (prominence_min != NULL) {
            bool* keep = _select_by_property_double(prominences, current_count, prominence_min, NULL);
            if (!keep) {
                free(all_peaks);
                free(left_edges);
                free(right_edges);
                free(current_peaks);
                free(prominences);
                free(left_bases);
                free(right_bases);
                free_peak_properties(properties);
                return 0;
            }
            
            // Update peak list
            size_t new_count = 0;
            for (size_t i = 0; i < current_count; i++) {
                if (keep[i]) {
                    current_peaks[new_count] = current_peaks[i];
                    prominences[new_count] = prominences[i];
                    left_bases[new_count] = left_bases[i];
                    right_bases[new_count] = right_bases[i];
                    if (properties->peak_heights) {
                        properties->peak_heights[new_count] = properties->peak_heights[i];
                    }
                    if (properties->plateau_sizes) {
                        properties->plateau_sizes[new_count] = properties->plateau_sizes[i];
                        properties->left_edges[new_count] = properties->left_edges[i];
                        properties->right_edges[new_count] = properties->right_edges[i];
                    }
                    new_count++;
                }
            }
            
            current_count = new_count;
            free(keep);
        }
        
        // Save prominence properties
        if (current_count > 0) {
            properties->prominences = (double*)malloc(current_count * sizeof(double));
            properties->left_bases = (size_t*)malloc(current_count * sizeof(size_t));
            properties->right_bases = (size_t*)malloc(current_count * sizeof(size_t));
            
            if (!properties->prominences || !properties->left_bases || !properties->right_bases) {
                free(all_peaks);
                free(left_edges);
                free(right_edges);
                free(current_peaks);
                free(prominences);
                free(left_bases);
                free(right_bases);
                free_peak_properties(properties);
                return 0;
            }
            
            memcpy(properties->prominences, prominences, current_count * sizeof(double));
            memcpy(properties->left_bases, left_bases, current_count * sizeof(size_t));
            memcpy(properties->right_bases, right_bases, current_count * sizeof(size_t));
        }
        
        // Evaluate width condition
        if (width_min != NULL) {
            // Calculate width
            double* widths = (double*)malloc(current_count * sizeof(double));
            double* width_heights = (double*)malloc(current_count * sizeof(double));
            double* left_ips = (double*)malloc(current_count * sizeof(double));
            double* right_ips = (double*)malloc(current_count * sizeof(double));
            
            if (!widths || !width_heights || !left_ips || !right_ips) {
                free(all_peaks);
                free(left_edges);
                free(right_edges);
                free(current_peaks);
                free(prominences);
                free(left_bases);
                free(right_bases);
                free(widths);
                free(width_heights);
                free(left_ips);
                free(right_ips);
                free_peak_properties(properties);
                return 0;
            }
            
            peak_widths(x, n, current_peaks, current_count, rel_height, prominences, 
                       left_bases, right_bases, widths, width_heights, left_ips, right_ips);
            
            // Filter peaks that meet width condition
            bool* keep = _select_by_property_double(widths, current_count, width_min, NULL);
            if (!keep) {
                free(all_peaks);
                free(left_edges);
                free(right_edges);
                free(current_peaks);
                free(prominences);
                free(left_bases);
                free(right_bases);
                free(widths);
                free(width_heights);
                free(left_ips);
                free(right_ips);
                free_peak_properties(properties);
                return 0;
            }
            
            // Update peak list
            size_t new_count = 0;
            for (size_t i = 0; i < current_count; i++) {
                if (keep[i]) {
                    current_peaks[new_count] = current_peaks[i];
                    widths[new_count] = widths[i];
                    width_heights[new_count] = width_heights[i];
                    left_ips[new_count] = left_ips[i];
                    right_ips[new_count] = right_ips[i];
                    if (properties->prominences) {
                        properties->prominences[new_count] = properties->prominences[i];
                        properties->left_bases[new_count] = properties->left_bases[i];
                        properties->right_bases[new_count] = properties->right_bases[i];
                    }
                    if (properties->peak_heights) {
                        properties->peak_heights[new_count] = properties->peak_heights[i];
                    }
                    if (properties->plateau_sizes) {
                        properties->plateau_sizes[new_count] = properties->plateau_sizes[i];
                        properties->left_edges[new_count] = properties->left_edges[i];
                        properties->right_edges[new_count] = properties->right_edges[i];
                    }
                    new_count++;
                }
            }
            
            current_count = new_count;
            
            // Save width properties
            if (current_count > 0) {
                properties->widths = (double*)malloc(current_count * sizeof(double));
                properties->width_heights = (double*)malloc(current_count * sizeof(double));
                properties->left_ips = (double*)malloc(current_count * sizeof(double));
                properties->right_ips = (double*)malloc(current_count * sizeof(double));
                
                if (!properties->widths || !properties->width_heights || 
                    !properties->left_ips || !properties->right_ips) {
                    free(all_peaks);
                    free(left_edges);
                    free(right_edges);
                    free(current_peaks);
                    free(prominences);
                    free(left_bases);
                    free(right_bases);
                    free(widths);
                    free(width_heights);
                    free(left_ips);
                    free(right_ips);
                    free(keep);
                    free_peak_properties(properties);
                    return 0;
                }
                
                memcpy(properties->widths, widths, current_count * sizeof(double));
                memcpy(properties->width_heights, width_heights, current_count * sizeof(double));
                memcpy(properties->left_ips, left_ips, current_count * sizeof(double));
                memcpy(properties->right_ips, right_ips, current_count * sizeof(double));
            }
            
            free(widths);
            free(width_heights);
            free(left_ips);
            free(right_ips);
            free(keep);
        }
        
        free(prominences);
        free(left_bases);
        free(right_bases);
    }
    
    // Return final peak list
    properties->count = current_count;
    *peaks = (size_t*)malloc(current_count * sizeof(size_t));
    if (!*peaks) {
        free(all_peaks);
        free(left_edges);
        free(right_edges);
        free(current_peaks);
        free_peak_properties(properties);
        return 0;
    }
    
    memcpy(*peaks, current_peaks, current_count * sizeof(size_t));
    
    // Free temporary arrays
    free(all_peaks);
    free(left_edges);
    free(right_edges);
    free(current_peaks);
    
    return current_count;
} 