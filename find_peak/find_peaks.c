#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <float.h>

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

// Initialize PeakResults structure
PeakResults* initialize_peak_results(int initial_capacity) {
    PeakResults* results = (PeakResults*)malloc(sizeof(PeakResults));
    if (!results) {
        return NULL;
    }
    
    results->capacity = initial_capacity;
    results->size = 0;
    
    results->indices = (int*)malloc(initial_capacity * sizeof(int));
    
    // Initialize other properties to NULL, to be allocated later as needed
    results->prominences = NULL;
    results->left_bases = NULL;
    results->right_bases = NULL;
    
    results->widths = NULL;
    results->width_heights = NULL;
    results->left_ips = NULL;
    results->right_ips = NULL;
    
    results->plateau_sizes = NULL;
    results->left_edges = NULL;
    results->right_edges = NULL;
    
    if (!results->indices) {
        free(results);
        return NULL;
    }
    
    return results;
}

// Free PeakResults structure and all allocated resources
void free_peak_results(PeakResults* results) {
    if (!results) {
        return;
    }
    
    free(results->indices);
    
    if (results->prominences) free(results->prominences);
    if (results->left_bases) free(results->left_bases);
    if (results->right_bases) free(results->right_bases);
    
    if (results->widths) free(results->widths);
    if (results->width_heights) free(results->width_heights);
    if (results->left_ips) free(results->left_ips);
    if (results->right_ips) free(results->right_ips);
    
    if (results->plateau_sizes) free(results->plateau_sizes);
    if (results->left_edges) free(results->left_edges);
    if (results->right_edges) free(results->right_edges);
    
    free(results);
}

// Increase results capacity
bool expand_peak_results(PeakResults* results) {
    int new_capacity = results->capacity * 2;
    int* new_indices = (int*)realloc(results->indices, new_capacity * sizeof(int));
    
    if (!new_indices) {
        return false;
    }
    
    results->indices = new_indices;
    results->capacity = new_capacity;
    
    return true;
}

// Add a peak
bool add_peak(PeakResults* results, int index) {
    if (results->size >= results->capacity) {
        if (!expand_peak_results(results)) {
            return false;
        }
    }
    
    results->indices[results->size++] = index;
    return true;
}

// Identify all local maxima
int _identify_local_maxima(const uint32_t* x, int size, PeakResults* results) {
    // Initialize peak results
    for (int i = 1; i < size - 1; i++) {
        // Handle flat peaks
        if (x[i] > x[i-1] && x[i] >= x[i+1]) {
            int j = i;
            while (j < size - 1 && x[j] == x[j+1]) {
                j++;
            }
            
            if (j == i) {
                // Single point peak
                add_peak(results, i);
            } else {
                // Flat peak, take the middle point
                add_peak(results, i + (j - i) / 2);
                i = j;
            }
        }
    }
    
    return results->size;
}

// Filter peaks by distance
void _filter_by_distance(PeakResults* results, const uint32_t* x, int distance) {
    if (results->size <= 1 || distance <= 1) {
        return;
    }
    
    // Sort by peak height first
    int* sorted_indices = (int*)malloc(results->size * sizeof(int));
    int* keep = (int*)calloc(results->size, sizeof(int));
    
    for (int i = 0; i < results->size; i++) {
        sorted_indices[i] = i;
    }
    
    // Sort peak indices by height from high to low
    for (int i = 0; i < results->size - 1; i++) {
        for (int j = 0; j < results->size - i - 1; j++) {
            if (x[results->indices[sorted_indices[j]]] < x[results->indices[sorted_indices[j+1]]]) {
                int temp = sorted_indices[j];
                sorted_indices[j] = sorted_indices[j+1];
                sorted_indices[j+1] = temp;
            }
        }
    }
    
    // Mark peaks to keep
    for (int i = 0; i < results->size; i++) {
        int idx = sorted_indices[i];
        
        if (!keep[idx]) {
            keep[idx] = 1;
            
            // Remove peaks that are too close
            for (int j = 0; j < results->size; j++) {
                if (!keep[j] && abs(results->indices[idx] - results->indices[j]) < distance) {
                    keep[j] = -1;  // Mark for removal
                }
            }
        }
    }
    
    // Reorganize the kept peaks
    int new_size = 0;
    for (int i = 0; i < results->size; i++) {
        if (keep[i] == 1) {
            results->indices[new_size++] = results->indices[i];
        }
    }
    
    results->size = new_size;
    
    // Sort by index
    for (int i = 0; i < results->size - 1; i++) {
        for (int j = 0; j < results->size - i - 1; j++) {
            if (results->indices[j] > results->indices[j+1]) {
                int temp = results->indices[j];
                results->indices[j] = results->indices[j+1];
                results->indices[j+1] = temp;
            }
        }
    }
    
    free(sorted_indices);
    free(keep);
}

// Calculate peak prominences
void _calculate_prominences(PeakResults* results, const uint32_t* x, int size, int wlen) {
    if (results->size == 0) {
        return;
    }
    
    results->prominences = (uint32_t*)malloc(results->size * sizeof(uint32_t));
    results->left_bases = (int*)malloc(results->size * sizeof(int));
    results->right_bases = (int*)malloc(results->size * sizeof(int));
    
    if (!results->prominences || !results->left_bases || !results->right_bases) {
        return;
    }
    
    for (int i = 0; i < results->size; i++) {
        int peak_idx = results->indices[i];
        uint32_t peak_height = x[peak_idx];
        
        // Determine search range
        int left_idx = (wlen > 0) ? (peak_idx - wlen / 2 > 0 ? peak_idx - wlen / 2 : 0) : 0;
        int right_idx = (wlen > 0) ? (peak_idx + wlen / 2 < size - 1 ? peak_idx + wlen / 2 : size - 1) : size - 1;
        
        // Find left and right minimum
        uint32_t left_min = peak_height;
        int left_min_idx = peak_idx;
        for (int j = peak_idx; j >= left_idx; j--) {
            if (x[j] < left_min) {
                left_min = x[j];
                left_min_idx = j;
            }
            
            // Stop if we encounter a higher peak
            if (j < peak_idx && x[j] > peak_height) {
                break;
            }
        }
        
        uint32_t right_min = peak_height;
        int right_min_idx = peak_idx;
        for (int j = peak_idx; j <= right_idx; j++) {
            if (x[j] < right_min) {
                right_min = x[j];
                right_min_idx = j;
            }
            
            // Stop if we encounter a higher peak
            if (j > peak_idx && x[j] > peak_height) {
                break;
            }
        }
        
        // Use the higher of the two minimums as reference level
        uint32_t reference_level = (left_min > right_min) ? left_min : right_min;
        
        // Calculate prominence
        results->prominences[i] = peak_height - reference_level;
        results->left_bases[i] = left_min_idx;
        results->right_bases[i] = right_min_idx;
    }
}

// Filter peaks by prominence
void _filter_by_prominence(PeakResults* results, uint32_t min_prominence, uint32_t max_prominence) {
    if (!results->prominences) {
        return;
    }
    
    int new_size = 0;
    for (int i = 0; i < results->size; i++) {
        if (results->prominences[i] >= min_prominence && 
            (max_prominence == 0 || results->prominences[i] <= max_prominence)) {
            
            results->indices[new_size] = results->indices[i];
            results->prominences[new_size] = results->prominences[i];
            results->left_bases[new_size] = results->left_bases[i];
            results->right_bases[new_size] = results->right_bases[i];
            new_size++;
        }
    }
    
    results->size = new_size;
}

// Calculate peak widths
void _calculate_widths(PeakResults* results, const uint32_t* x, int size, uint32_t rel_height) {
    if (results->size == 0 || !results->prominences) {
        return;
    }
    
    results->widths = (uint32_t*)malloc(results->size * sizeof(uint32_t));
    results->width_heights = (uint32_t*)malloc(results->size * sizeof(uint32_t));
    results->left_ips = (uint32_t*)malloc(results->size * sizeof(uint32_t));
    results->right_ips = (uint32_t*)malloc(results->size * sizeof(uint32_t));
    
    if (!results->widths || !results->width_heights || !results->left_ips || !results->right_ips) {
        return;
    }
    
    for (int i = 0; i < results->size; i++) {
        int peak_idx = results->indices[i];
        uint32_t peak_height = x[peak_idx];
        
        // Calculate reference height
        uint32_t height = peak_height - (results->prominences[i] * rel_height) / 100;
        results->width_heights[i] = height;
        
        // Find left intersection point
        int left_idx = peak_idx;
        while (left_idx > 0 && x[left_idx] >= height) {
            left_idx--;
        }
        
        // Linear interpolation for precise intersection point
        uint32_t left_ip;
        if (left_idx == peak_idx) {
            left_ip = left_idx;
        } else if (x[left_idx] == height) {
            left_ip = left_idx;
        } else if (x[left_idx + 1] != x[left_idx]) {
            left_ip = left_idx + (height - x[left_idx]) * 100 / (x[left_idx + 1] - x[left_idx]);
        } else {
            left_ip = left_idx;
        }
        results->left_ips[i] = left_ip;
        
        // Find right intersection point
        int right_idx = peak_idx;
        while (right_idx < size - 1 && x[right_idx] >= height) {
            right_idx++;
        }
        
        // Linear interpolation for precise intersection point
        uint32_t right_ip;
        if (right_idx == peak_idx) {
            right_ip = right_idx;
        } else if (x[right_idx] == height) {
            right_ip = right_idx;
        } else if (x[right_idx] != x[right_idx - 1]) {
            right_ip = right_idx - 1 + (height - x[right_idx - 1]) * 100 / (x[right_idx] - x[right_idx - 1]);
        } else {
            right_ip = right_idx;
        }
        results->right_ips[i] = right_ip;
        
        // Calculate width
        results->widths[i] = right_ip - left_ip;
    }
}

// Filter peaks by width
void _filter_by_width(PeakResults* results, uint32_t min_width, uint32_t max_width) {
    if (!results->widths) {
        return;
    }
    
    int new_size = 0;
    for (int i = 0; i < results->size; i++) {
        if (results->widths[i] >= min_width && 
            (max_width == 0 || results->widths[i] <= max_width)) {
            
            // Copy all properties
            results->indices[new_size] = results->indices[i];
            
            if (results->prominences) {
                results->prominences[new_size] = results->prominences[i];
                results->left_bases[new_size] = results->left_bases[i];
                results->right_bases[new_size] = results->right_bases[i];
            }
            
            results->widths[new_size] = results->widths[i];
            results->width_heights[new_size] = results->width_heights[i];
            results->left_ips[new_size] = results->left_ips[i];
            results->right_ips[new_size] = results->right_ips[i];
            
            if (results->plateau_sizes) {
                results->plateau_sizes[new_size] = results->plateau_sizes[i];
                results->left_edges[new_size] = results->left_edges[i];
                results->right_edges[new_size] = results->right_edges[i];
            }
            
            new_size++;
        }
    }
    
    results->size = new_size;
}

// Calculate plateau sizes
void _calculate_plateau_sizes(PeakResults* results, const uint32_t* x, int size) {
    if (results->size == 0) {
        return;
    }
    
    results->plateau_sizes = (int*)malloc(results->size * sizeof(int));
    results->left_edges = (int*)malloc(results->size * sizeof(int));
    results->right_edges = (int*)malloc(results->size * sizeof(int));
    
    if (!results->plateau_sizes || !results->left_edges || !results->right_edges) {
        return;
    }
    
    for (int i = 0; i < results->size; i++) {
        int peak_idx = results->indices[i];
        uint32_t peak_height = x[peak_idx];
        
        // Find left plateau edge
        int left_edge = peak_idx;
        while (left_edge > 0 && x[left_edge - 1] == peak_height) {
            left_edge--;
        }
        results->left_edges[i] = left_edge;
        
        // Find right plateau edge
        int right_edge = peak_idx;
        while (right_edge < size - 1 && x[right_edge + 1] == peak_height) {
            right_edge++;
        }
        results->right_edges[i] = right_edge;
        
        // Calculate plateau size
        results->plateau_sizes[i] = right_edge - left_edge + 1;
    }
}

// Filter peaks by plateau size
void _filter_by_plateau_size(PeakResults* results, int min_plateau_size, int max_plateau_size) {
    if (!results->plateau_sizes) {
        return;
    }
    
    int new_size = 0;
    for (int i = 0; i < results->size; i++) {
        if (results->plateau_sizes[i] >= min_plateau_size && 
            (max_plateau_size <= 0 || results->plateau_sizes[i] <= max_plateau_size)) {
            
            // Copy all properties
            results->indices[new_size] = results->indices[i];
            
            if (results->prominences) {
                results->prominences[new_size] = results->prominences[i];
                results->left_bases[new_size] = results->left_bases[i];
                results->right_bases[new_size] = results->right_bases[i];
            }
            
            if (results->widths) {
                results->widths[new_size] = results->widths[i];
                results->width_heights[new_size] = results->width_heights[i];
                results->left_ips[new_size] = results->left_ips[i];
                results->right_ips[new_size] = results->right_ips[i];
            }
            
            results->plateau_sizes[new_size] = results->plateau_sizes[i];
            results->left_edges[new_size] = results->left_edges[i];
            results->right_edges[new_size] = results->right_edges[i];
            
            new_size++;
        }
    }
    
    results->size = new_size;
}

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
                        int wlen, uint32_t rel_height) {
    
    if (!x || size <= 2) {
        return NULL;
    }
    
    // Initialize peak results structure
    PeakResults* results = initialize_peak_results(size / 10 > 10 ? size / 10 : 10);
    if (!results) {
        return NULL;
    }
    
    // Find all local maxima
    _identify_local_maxima(x, size, results);
    
    // Apply plateau size filtering
    if (min_plateau_size > 0 || max_plateau_size > 0) {
        _calculate_plateau_sizes(results, x, size);
        _filter_by_plateau_size(results, min_plateau_size, max_plateau_size);
    }
    
    // Apply distance filtering
    if (distance > 1) {
        _filter_by_distance(results, x, distance);
    }
    
    // Apply prominence filtering
    if (min_prominence > 0 || max_prominence > 0) {
        _calculate_prominences(results, x, size, wlen);
        _filter_by_prominence(results, min_prominence, max_prominence);
    }
    
    // Apply width filtering
    if (min_width > 0 || max_width > 0) {
        if (!results->prominences) {
            _calculate_prominences(results, x, size, wlen);
        }
        _calculate_widths(results, x, size, rel_height);
        _filter_by_width(results, min_width, max_width);
    }
    
    return results;
}

// Simple example program that can be called separately
void run_example() {
    // Create a simple signal
    int size = 100;
    uint32_t* signal = (uint32_t*)malloc(size * sizeof(uint32_t));
    
    // Generate a signal with several peaks
    for (int i = 0; i < size; i++) {
        signal[i] = 0;
        
        // Add a few Gaussian curves as peaks
        signal[i] += (uint32_t)(5 * exp(-0.1 * (i - 20) * (i - 20)));  // Peak at i=20
        signal[i] += (uint32_t)(3 * exp(-0.1 * (i - 50) * (i - 50)));  // Peak at i=50
        signal[i] += (uint32_t)(7 * exp(-0.1 * (i - 80) * (i - 80)));  // Peak at i=80
    }
    
    // Call find_peaks function
    PeakResults* results = find_peaks(signal, size, 
                                      10,    // distance
                                      2, 0,  // min_prominence, max_prominence
                                      3, 0,  // min_width, max_width
                                      0, 0,  // min_plateau_size, max_plateau_size
                                      0,     // wlen
                                      50);   // rel_height (50 = 0.5 * 100)
    
    if (results) {
        printf("Found %d peaks:\n", results->size);
        for (int i = 0; i < results->size; i++) {
            printf("Peak #%d: index=%d, height=%u", i+1, results->indices[i], signal[results->indices[i]]);
            
            if (results->prominences) {
                printf(", prominence=%u", results->prominences[i]);
            }
            
            if (results->widths) {
                printf(", width=%u", results->widths[i]);
            }
            
            if (results->plateau_sizes) {
                printf(", plateau_size=%d", results->plateau_sizes[i]);
            }
            
            printf("\n");
        }
        
        // Free results and signal memory
        free_peak_results(results);
    } else {
        printf("No peaks found or error occurred during processing\n");
    }
    
    free(signal);
}

/* 
// Main function removed to avoid conflicts with test_peaks.c
int main() {
    run_example();
    return 0;
}
*/ 