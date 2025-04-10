#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#define MAX_SIGNAL_LENGTH 1000
#define MAX_REGIONS 100
#define MAX_FILENAME_LENGTH 256

typedef struct {
    uint32_t start;
    uint32_t end;
    uint32_t values[MAX_SIGNAL_LENGTH];
    int values_count;
} TargetRegion;

typedef struct {
    uint32_t start;
    uint32_t end;
    char trend_type[20]; // "Upward Trend" or "Downward Trend"
} TrendRegion;

/**
 * Detect oscillation regions in a signal
 * 
 * Parameters:
 * - signal: Array of signal values
 * - signal_length: Length of the signal array
 * - min_length: Minimum length of the interval
 * - debug: Whether to output debug information
 * - signal_name: Signal name (for debugging)
 * - equal_tolerance: Tolerance for judging value equality
 * - target_regions: Output array to store detected target regions
 * - trend_regions: Output array to store detected trend regions
 * - trend_count: Pointer to store the number of trend regions detected
 * 
 * Returns:
 * - Number of target regions detected
 */
int detect_oscillation_regions(
    uint32_t* signal,
    int signal_length,
    int min_length,
    bool debug,
    const char* signal_name,
    uint32_t equal_tolerance,
    TargetRegion* target_regions,
    TrendRegion* trend_regions,
    int* trend_count_ptr
) {
    if (signal_length == 0) {
        return 0;
    }
    
    // Find valid indices (non-zero values)
    int valid_indices[MAX_SIGNAL_LENGTH];
    uint32_t valid_signal[MAX_SIGNAL_LENGTH];
    int valid_count = 0;
    
    for (int i = 0; i < signal_length; i++) {
        if (signal[i] > 0) {
            valid_indices[valid_count] = i;
            valid_signal[valid_count] = signal[i];
            valid_count++;
        }
    }
    
    if (valid_count == 0) {
        return 0;
    }
    
    if (debug && signal_name != NULL) {
        printf("\n===== %s Detailed Debug Info =====\n", signal_name);
        printf("Valid data points: %d\n", valid_count);
        
        // Count unique values (simple implementation)
        int unique_count = 1;
        for (int i = 1; i < valid_count; i++) {
            bool is_unique = true;
            for (int j = 0; j < i; j++) {
                if (valid_signal[i] == valid_signal[j]) {
                    is_unique = false;
                    break;
                }
            }
            if (is_unique) {
                unique_count++;
            }
        }
        printf("Unique values: %d\n", unique_count);
    }
    
    // Calculate direction changes between points
    char directions[MAX_SIGNAL_LENGTH][5]; // "up", "down", or "flat"
    for (int i = 1; i < valid_count; i++) {
        if (valid_signal[i] > valid_signal[i-1] + equal_tolerance) {
            strcpy(directions[i-1], "up");
        } else if (valid_signal[i] < valid_signal[i-1] - equal_tolerance) {
            strcpy(directions[i-1], "down");
        } else {
            strcpy(directions[i-1], "flat");
        }
    }
    
    // Identify upward and downward trend intervals
    int trend_count = 0;
    char current_direction[5] = "";
    int start_idx = 0;
    
    for (int i = 0; i < valid_count - 1; i++) {
        // To align indices with original valid_indices, add 1
        int actual_idx = i + 1;
        
        // Start a new trend interval
        if (strlen(current_direction) == 0 && 
            (strcmp(directions[i], "up") == 0 || strcmp(directions[i], "down") == 0)) {
            strcpy(current_direction, directions[i]);
            start_idx = i;
        }
        // Direction change
        else if (strlen(current_direction) > 0 && 
                 strcmp(directions[i], "flat") != 0 && 
                 strcmp(directions[i], current_direction) != 0) {
            // Check interval length
            if (actual_idx - start_idx >= min_length) {
                uint32_t region_start = valid_indices[start_idx];
                uint32_t region_end = valid_indices[actual_idx - 1];
                
                if (trend_count < MAX_REGIONS) {
                    trend_regions[trend_count].start = region_start;
                    trend_regions[trend_count].end = region_end;
                    
                    if (strcmp(current_direction, "up") == 0) {
                        strcpy(trend_regions[trend_count].trend_type, "Upward Trend");
                    } else {
                        strcpy(trend_regions[trend_count].trend_type, "Downward Trend");
                    }
                    
                    if (debug) {
                        printf("Detected %s: from time point %u(a) to %u(b), length: %u\n",
                               trend_regions[trend_count].trend_type,
                               region_start, region_end, region_end - region_start + 1);
                    }
                    
                    trend_count++;
                }
            }
            
            // Start new trend interval
            strcpy(current_direction, directions[i]);
            start_idx = i;
        }
    }
    
    // Process the last trend interval
    if (strlen(current_direction) > 0 && valid_count - start_idx >= min_length) {
        uint32_t region_start = valid_indices[start_idx];
        uint32_t region_end = valid_indices[valid_count - 1];
        
        if (trend_count < MAX_REGIONS) {
            trend_regions[trend_count].start = region_start;
            trend_regions[trend_count].end = region_end;
            
            if (strcmp(current_direction, "up") == 0) {
                strcpy(trend_regions[trend_count].trend_type, "Upward Trend");
            } else {
                strcpy(trend_regions[trend_count].trend_type, "Downward Trend");
            }
            
            if (debug) {
                printf("Detected %s: from time point %u(a) to %u(b), length: %u\n",
                       trend_regions[trend_count].trend_type,
                       region_start, region_end, region_end - region_start + 1);
            }
            
            trend_count++;
        }
    }
    
    // Store trend count for later use
    *trend_count_ptr = trend_count;
    
    // Find target intervals (intervals not belonging to trends) based on trend intervals
    int target_count = 0;
    
    // Sort trend intervals by start time (simple bubble sort)
    for (int i = 0; i < trend_count - 1; i++) {
        for (int j = 0; j < trend_count - i - 1; j++) {
            if (trend_regions[j].start > trend_regions[j + 1].start) {
                TrendRegion temp = trend_regions[j];
                trend_regions[j] = trend_regions[j + 1];
                trend_regions[j + 1] = temp;
            }
        }
    }
    
    // If there are trend intervals, check gaps between trend intervals
    if (trend_count > 0) {
        // Check gap before the first trend interval
        if (trend_regions[0].start > (uint32_t)valid_indices[0] && 
            trend_regions[0].start - (uint32_t)valid_indices[0] + 1 >= (uint32_t)min_length) {
            
            uint32_t target_start = (uint32_t)valid_indices[0];
            uint32_t target_end = trend_regions[0].start - 1;
            
            if (debug) {
                printf("Detected target region: from time point %u(s) to %u(t), length: %u\n",
                       target_start, target_end, target_end - target_start + 1);
            }
            
            target_regions[target_count].start = target_start;
            target_regions[target_count].end = target_end;
            target_regions[target_count].values_count = 0;
            
            // Copy non-zero values
            for (uint32_t i = target_start; i <= target_end; i++) {
                if (signal[i] > 0) {
                    target_regions[target_count].values[target_regions[target_count].values_count++] = signal[i];
                }
            }
            
            target_count++;
        }
        
        // Check gaps between trend intervals
        for (int i = 0; i < trend_count - 1; i++) {
            uint32_t current_end = trend_regions[i].end;
            uint32_t next_start = trend_regions[i+1].start;
            
            if (next_start > current_end + 1 && next_start - current_end - 1 >= (uint32_t)min_length) {
                uint32_t target_start = current_end + 1;
                uint32_t target_end = next_start - 1;
                
                if (debug) {
                    printf("Detected target region: from time point %u(s) to %u(t), length: %u\n",
                           target_start, target_end, target_end - target_start + 1);
                }
                
                target_regions[target_count].start = target_start;
                target_regions[target_count].end = target_end;
                target_regions[target_count].values_count = 0;
                
                // Copy non-zero values
                for (uint32_t i = target_start; i <= target_end; i++) {
                    if (signal[i] > 0) {
                        target_regions[target_count].values[target_regions[target_count].values_count++] = signal[i];
                    }
                }
                
                target_count++;
            }
        }
        
        // Check gap after the last trend interval
        if (trend_regions[trend_count-1].end < (uint32_t)valid_indices[valid_count-1] && 
            (uint32_t)valid_indices[valid_count-1] - trend_regions[trend_count-1].end >= (uint32_t)min_length) {
            
            uint32_t target_start = trend_regions[trend_count-1].end + 1;
            uint32_t target_end = (uint32_t)valid_indices[valid_count-1];
            
            if (debug) {
                printf("Detected target region: from time point %u(s) to %u(t), length: %u\n",
                       target_start, target_end, target_end - target_start + 1);
            }
            
            target_regions[target_count].start = target_start;
            target_regions[target_count].end = target_end;
            target_regions[target_count].values_count = 0;
            
            // Copy non-zero values
            for (uint32_t i = target_start; i <= target_end; i++) {
                if (signal[i] > 0) {
                    target_regions[target_count].values[target_regions[target_count].values_count++] = signal[i];
                }
            }
            
            target_count++;
        }
    }
    // If there are no trend intervals, check if the entire signal can be considered as a target interval
    else if (valid_count >= min_length) {
        uint32_t target_start = valid_indices[0];
        uint32_t target_end = valid_indices[valid_count-1];
        
        if (debug) {
            printf("Detected target region (no trend regions): from time point %u(s) to %u(t), length: %u\n",
                   target_start, target_end, target_end - target_start + 1);
        }
        
        target_regions[target_count].start = target_start;
        target_regions[target_count].end = target_end;
        target_regions[target_count].values_count = 0;
        
        // Copy valid values
        for (int i = 0; i < valid_count; i++) {
            target_regions[target_count].values[i] = valid_signal[i];
        }
        target_regions[target_count].values_count = valid_count;
        
        target_count++;
    }
    
    if (debug && signal_name != NULL) {
        printf("\nNumber of target regions detected: %d\n", target_count);
        for (int i = 0; i < target_count; i++) {
            printf("  Region %d: time points %u(s)-%u(t) (length: %u)\n",
                   i+1, target_regions[i].start, target_regions[i].end,
                   target_regions[i].end - target_regions[i].start + 1);
        }
    }
    
    return target_count;
}

/**
 * Load signal data from a CSV file
 * 
 * Parameters:
 * - filename: Name of the CSV file
 * - x_values: Array to store x values
 * - signals: 2D array to store signal values (signals[signal_index][time_index])
 * - max_signals: Maximum number of signals to read
 * - signal_names: Array to store signal names
 * 
 * Returns:
 * - Number of data points (rows) read
 */
int load_signal_data(
    const char* filename,
    uint32_t* x_values,
    uint32_t signals[][MAX_SIGNAL_LENGTH],
    int max_signals,
    char signal_names[][MAX_FILENAME_LENGTH],
    int* signal_count
) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Data file not found. Please ensure %s exists\n", filename);
        return 0;
    }
    
    printf("Using data: %s\n", filename);
    
    char line[1024];
    int row_count = 0;
    *signal_count = 0;
    
    // Read header to get signal names
    if (fgets(line, sizeof(line), file)) {
        char* token = strtok(line, ",");
        int col = 0;
        
        // Skip "x" column
        if (token != NULL) {
            token = strtok(NULL, ",");
            col++;
        }
        
        // Read signal names
        while (token != NULL && col <= max_signals) {
            // Remove newline if present
            size_t len = strlen(token);
            if (len > 0 && (token[len-1] == '\n' || token[len-1] == '\r')) {
                token[len-1] = '\0';
            }
            
            strncpy(signal_names[col-1], token, MAX_FILENAME_LENGTH-1);
            signal_names[col-1][MAX_FILENAME_LENGTH-1] = '\0';
            
            token = strtok(NULL, ",");
            col++;
        }
        
        *signal_count = col - 1;
    }
    
    // Read data rows
    while (fgets(line, sizeof(line), file) && row_count < MAX_SIGNAL_LENGTH) {
        char* token = strtok(line, ",");
        int col = 0;
        
        while (token != NULL && col <= *signal_count) {
            uint32_t value = (uint32_t)atoi(token);
            
            if (col == 0) {
                x_values[row_count] = value;
            } else {
                signals[col-1][row_count] = value;
            }
            
            token = strtok(NULL, ",");
            col++;
        }
        
        row_count++;
    }
    
    fclose(file);
    return row_count;
}

/**
 * Generate a Gnuplot script to visualize the data and detected regions
 * 
 * Parameters:
 * - filename: Name of the output Gnuplot script file
 * - data_file: Name of the CSV data file
 * - x_values: Array of x values
 * - signals: 2D array of signal values
 * - signal_count: Number of signals
 * - data_count: Number of data points
 * - signal_names: Array of signal names
 * - target_regions: Array of target regions for each signal
 * - target_counts: Array of target region counts for each signal
 * - trend_regions: Array of trend regions for each signal
 * - trend_counts: Array of trend region counts for each signal
 * - min_length: Minimum length of intervals
 * 
 * Returns:
 * - true if successful, false otherwise
 */
bool generate_gnuplot_script(
    const char* filename,
    const char* data_file,
    uint32_t* x_values,
    uint32_t signals[][MAX_SIGNAL_LENGTH],
    int signal_count,
    int data_count,
    char signal_names[][MAX_FILENAME_LENGTH],
    TargetRegion target_regions[][MAX_REGIONS],
    int* target_counts,
    TrendRegion trend_regions[][MAX_REGIONS],
    int* trend_counts,
    int min_length
) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Failed to create Gnuplot script file: %s\n", filename);
        return false;
    }
    
    // Write Gnuplot script header
    fprintf(file, "# Gnuplot script for oscillation visualization\n");
    fprintf(file, "set terminal png size 1600,1200 enhanced font 'Arial,12'\n");
    fprintf(file, "set output 'signal_visualization.png'\n");
    fprintf(file, "set title 'Signal Visualization - Target Region Markers (s,t=Start/End points, Length≥%d)' font 'Arial,16'\n", min_length);
    fprintf(file, "set xlabel 'Time Points' font 'Arial,12'\n");
    fprintf(file, "set ylabel 'Signal Values' font 'Arial,12'\n");
    fprintf(file, "set grid\n");
    fprintf(file, "set key top right\n");
    
    // Calculate y-range
    uint32_t max_value = 0;
    for (int i = 0; i < signal_count; i++) {
        for (int j = 0; j < data_count; j++) {
            if (signals[i][j] > max_value) {
                max_value = signals[i][j];
            }
        }
    }
    fprintf(file, "set yrange [0:%d]\n", (int)(max_value * 1.3));
    
    // Generate plot command
    fprintf(file, "plot \\\n");
    
    // Plot each signal
    const char* colors[] = {"red", "blue", "green", "purple", "orange", "brown", "black", "cyan", "magenta", "yellow"};
    
    // First, create a temporary data file for each signal
    for (int i = 0; i < signal_count; i++) {
        char temp_file[256];
        sprintf(temp_file, "signal_%d.tmp", i);
        
        FILE* temp = fopen(temp_file, "w");
        if (temp == NULL) {
            fprintf(stderr, "Failed to create temporary data file: %s\n", temp_file);
            fclose(file);
            return false;
        }
        
        // Write valid data points
        for (int j = 0; j < data_count; j++) {
            if (signals[i][j] > 0) {
                fprintf(temp, "%d %d\n", x_values[j], signals[i][j]);
            }
        }
        
        fclose(temp);
        
        // Add plot command for this signal
        int color_idx = i % 10;
        fprintf(file, "    'signal_%d.tmp' using 1:2 with linespoints lt 1 pt 7 ps 0.5 lc rgb '%s' title '%s'", 
                i, colors[color_idx], signal_names[i]);
        
        // Add vertical lines for target regions
        for (int j = 0; j < target_counts[i]; j++) {
            fprintf(file, ", \\\n    %d with lines lt 0 lw 2 lc rgb '%s' notitle", 
                    target_regions[i][j].start, colors[color_idx]);
            fprintf(file, ", \\\n    %d with lines lt 0 lw 2 lc rgb '%s' notitle", 
                    target_regions[i][j].end, colors[color_idx]);
            
            // Add 's' and 't' labels
            int y_pos = max_value + 3 + i * 2;
            fprintf(file, ", \\\n    '-' with labels point pt 7 ps 0.5 offset 0,0.5 tc rgb '%s' notitle", 
                    colors[color_idx]);
            fprintf(file, ", \\\n    '-' with labels point pt 7 ps 0.5 offset 0,0.5 tc rgb '%s' notitle", 
                    colors[color_idx]);
        }
        
        if (i < signal_count - 1) {
            fprintf(file, ", \\\n");
        } else {
            fprintf(file, "\n");
        }
    }
    
    // Add label data
    for (int i = 0; i < signal_count; i++) {
        for (int j = 0; j < target_counts[i]; j++) {
            int y_pos = max_value + 3 + i * 2;
            fprintf(file, "%d %d \"s\"\ne\n", target_regions[i][j].start, y_pos);
            fprintf(file, "%d %d \"t\"\ne\n", target_regions[i][j].end, y_pos);
        }
    }
    
    fclose(file);
    printf("Gnuplot script generated: %s\n", filename);
    return true;
}

/**
 * Run Gnuplot to generate the visualization
 * 
 * Parameters:
 * - script_file: Name of the Gnuplot script file
 * 
 * Returns:
 * - true if successful, false otherwise
 */
bool run_gnuplot(const char* script_file) {
    char command[512];
    
    // Check if Gnuplot is installed
    printf("Checking for Gnuplot installation...\n");
    FILE* test = popen("gnuplot --version", "r");
    if (test == NULL) {
        fprintf(stderr, "Failed to check for Gnuplot. Please make sure Gnuplot is installed.\n");
        return false;
    }
    
    char version[256];
    if (fgets(version, sizeof(version), test) == NULL) {
        fprintf(stderr, "Gnuplot not found. Please install Gnuplot to generate visualizations.\n");
        pclose(test);
        return false;
    }
    
    pclose(test);
    printf("Found %s", version);
    
    // Run Gnuplot
    printf("Running Gnuplot to generate visualization...\n");
    sprintf(command, "gnuplot \"%s\"", script_file);
    
    int result = system(command);
    if (result != 0) {
        fprintf(stderr, "Failed to run Gnuplot. Error code: %d\n", result);
        return false;
    }
    
    printf("Visualization generated successfully: signal_visualization.png\n");
    return true;
}

/**
 * Main function - entry point of the program
 */
int main() {
    const char* data_file = "refined_test_data.csv";
    uint32_t x_values[MAX_SIGNAL_LENGTH];
    uint32_t signals[10][MAX_SIGNAL_LENGTH]; // Support up to 10 signals
    char signal_names[10][MAX_FILENAME_LENGTH];
    int signal_count = 0;
    
    // Arrays to store detected regions for each signal
    TargetRegion target_regions[10][MAX_REGIONS];
    TrendRegion trend_regions[10][MAX_REGIONS];
    int target_counts[10] = {0};
    int trend_counts[10] = {0};
    
    // Define target interval detection parameters
    int min_length = 5;  // Minimum length of the target interval
    
    // Load signal data
    int data_count = load_signal_data(data_file, x_values, signals, 10, signal_names, &signal_count);
    
    if (data_count == 0) {
        return 1;
    }
    
    // Process each signal and detect oscillation regions
    for (int i = 0; i < signal_count; i++) {
        target_counts[i] = detect_oscillation_regions(
            signals[i],
            data_count,
            min_length,
            true,  // Enable debugging
            signal_names[i],
            1,     // Equal tolerance (converted from 0.25 to integer)
            target_regions[i],
            trend_regions[i],
            &trend_counts[i]
        );
        
        printf("Signal %s: detected %d oscillation regions\n", signal_names[i], target_counts[i]);
    }
    
    // Check if visualization should be skipped
    const char* skip_viz = getenv("SKIP_VISUALIZATION");
    if (skip_viz != NULL && strcmp(skip_viz, "1") == 0) {
        printf("\nSkipping visualization as requested.\n");
        return 0;
    }
    
    // Generate Gnuplot visualization
    const char* gnuplot_script = "plot_signals.gp";
    if (generate_gnuplot_script(
            gnuplot_script,
            data_file,
            x_values,
            signals,
            signal_count,
            data_count,
            signal_names,
            target_regions,
            target_counts,
            trend_regions,
            trend_counts,
            min_length)) {
        
        // Run Gnuplot to generate the visualization
        if (run_gnuplot(gnuplot_script)) {
            printf("\nVisualization completed successfully.\n");
        } else {
            printf("\nFailed to generate visualization. Please install Gnuplot and try again.\n");
            printf("You can manually run the generated script with: gnuplot %s\n", gnuplot_script);
        }
    }
    
    return 0;
} 