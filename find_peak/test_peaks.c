#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "find_peaks.h"

// Generate test signal
void generate_test_signal(uint32_t* signal, int size) {
    for (int i = 0; i < size; i++) {
        signal[i] = 0;
        
        // Add a few Gaussian curves as peaks
        signal[i] += (uint32_t)(5 * exp(-0.1 * (i - 20) * (i - 20)));  // Peak at i=20
        signal[i] += (uint32_t)(3 * exp(-0.1 * (i - 50) * (i - 50)));  // Peak at i=50
        signal[i] += (uint32_t)(7 * exp(-0.1 * (i - 80) * (i - 80)));  // Peak at i=80
        
        // Add some plateau peaks
        if (i >= 30 && i <= 35) signal[i] = 4;  // Plateau peak
        
        // Add some noise
        signal[i] += (uint32_t)(((double)rand() / RAND_MAX - 0.5) * 0.5);
    }
}

// Save signal to file for plotting
void save_signal_to_file(const uint32_t* signal, int size, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Cannot create file: %s\n", filename);
        return;
    }
    
    for (int i = 0; i < size; i++) {
        fprintf(file, "%d,%u\n", i, signal[i]);
    }
    
    fclose(file);
}

// Save peak results to file
void save_peaks_to_file(const uint32_t* signal, const PeakResults* results, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Cannot create file: %s\n", filename);
        return;
    }
    
    for (int i = 0; i < results->size; i++) {
        int idx = results->indices[i];
        fprintf(file, "%d,%u", idx, signal[idx]);
        
        if (results->prominences) {
            fprintf(file, ",%u", results->prominences[i]);
        }
        
        if (results->widths) {
            fprintf(file, ",%u", results->widths[i]);
        }
        
        if (results->plateau_sizes) {
            fprintf(file, ",%d", results->plateau_sizes[i]);
        }
        
        fprintf(file, "\n");
    }
    
    fclose(file);
}

int main() {
    // Set signal size
    int signal_size = 100;
    uint32_t* signal = (uint32_t*)malloc(signal_size * sizeof(uint32_t));
    
    // Generate test signal
    generate_test_signal(signal, signal_size);
    
    // Save signal to file
    save_signal_to_file(signal, signal_size, "signal.csv");
    
    printf("Testing distance parameter:\n");
    PeakResults* results1 = find_peaks(signal, signal_size, 
                                      15,   // distance = 15
                                      0, 0, // prominence not limited
                                      0, 0, // width not limited
                                      0, 0, // plateau_size not limited
                                      0, 50);
    
    if (results1) {
        printf("Found %d peaks using distance=15\n", results1->size);
        for (int i = 0; i < results1->size; i++) {
            printf("  Peak #%d: index=%d, height=%u\n", 
                   i+1, results1->indices[i], signal[results1->indices[i]]);
        }
        save_peaks_to_file(signal, results1, "peaks_distance.csv");
        free_peak_results(results1);
    }
    
    printf("\nTesting prominence parameter:\n");
    PeakResults* results2 = find_peaks(signal, signal_size, 
                                      0,    // distance not limited
                                      3, 0, // min_prominence = 3
                                      0, 0, // width not limited
                                      0, 0, // plateau_size not limited
                                      0, 50);
    
    if (results2) {
        printf("Found %d peaks using min_prominence=3\n", results2->size);
        for (int i = 0; i < results2->size; i++) {
            printf("  Peak #%d: index=%d, height=%u, prominence=%u\n", 
                   i+1, results2->indices[i], signal[results2->indices[i]], results2->prominences[i]);
        }
        save_peaks_to_file(signal, results2, "peaks_prominence.csv");
        free_peak_results(results2);
    }
    
    printf("\nTesting width parameter:\n");
    PeakResults* results3 = find_peaks(signal, signal_size, 
                                      0,    // distance not limited
                                      0, 0, // prominence not limited
                                      5, 0, // min_width = 5
                                      0, 0, // plateau_size not limited
                                      0, 50);
    
    if (results3) {
        printf("Found %d peaks using min_width=5\n", results3->size);
        for (int i = 0; i < results3->size; i++) {
            printf("  Peak #%d: index=%d, height=%u, width=%u\n", 
                   i+1, results3->indices[i], signal[results3->indices[i]], results3->widths[i]);
        }
        save_peaks_to_file(signal, results3, "peaks_width.csv");
        free_peak_results(results3);
    }
    
    printf("\nTesting plateau size parameter:\n");
    PeakResults* results4 = find_peaks(signal, signal_size, 
                                      0,    // distance not limited
                                      0, 0, // prominence not limited
                                      0, 0, // width not limited
                                      3, 0, // min_plateau_size = 3
                                      0, 50);
    
    if (results4) {
        printf("Found %d peaks using min_plateau_size=3\n", results4->size);
        for (int i = 0; i < results4->size; i++) {
            printf("  Peak #%d: index=%d, height=%u, plateau_size=%d\n", 
                   i+1, results4->indices[i], signal[results4->indices[i]], results4->plateau_sizes[i]);
        }
        save_peaks_to_file(signal, results4, "peaks_plateau.csv");
        free_peak_results(results4);
    }
    
    printf("\nTesting combined parameters:\n");
    PeakResults* results5 = find_peaks(signal, signal_size, 
                                      10,   // distance = 10
                                      2, 0, // min_prominence = 2
                                      3, 0, // min_width = 3
                                      0, 0, // plateau_size not limited
                                      0, 50);
    
    if (results5) {
        printf("Found %d peaks using combined parameters\n", results5->size);
        for (int i = 0; i < results5->size; i++) {
            printf("  Peak #%d: index=%d, height=%u, prominence=%u, width=%u\n", 
                   i+1, results5->indices[i], signal[results5->indices[i]], 
                   results5->prominences[i], results5->widths[i]);
        }
        save_peaks_to_file(signal, results5, "peaks_combined.csv");
        free_peak_results(results5);
    }
    
    // Free signal memory
    free(signal);
    
    return 0;
} 