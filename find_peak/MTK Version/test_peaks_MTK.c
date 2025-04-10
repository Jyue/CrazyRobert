#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "find_peaks_MTK.h"

// 當在公司環境中使用時，這兩個函數將被實際的MTK對應函數替代
// 這裡只是為了編譯能通過的暫時替代品
void* get_ctrl_buffer(unsigned int size) {
    // 在實際環境中，這將調用真正的get_ctrl_buffer函數
    printf("Called get_ctrl_buffer(%u)\n", size);
    return malloc(size);
}

void free_ctrl_buffer(void* buff_ptr) {
    // 在實際環境中，這將調用真正的free_ctrl_buffer函數
    printf("Called free_ctrl_buffer(0x%p)\n", buff_ptr);
    free(buff_ptr);
}

// 生成測試信號
uint32_t* generate_test_signal(int* signal_size) {
    int size = 100;
    *signal_size = size;
    
    uint32_t* signal = (uint32_t*)get_ctrl_buffer(size * sizeof(uint32_t));
    
    // 生成一個具有幾個峰值的信號
    for (int i = 0; i < size; i++) {
        signal[i] = 0;
        
        // 加入幾個高斯曲線作為峰值
        signal[i] += (uint32_t)(5 * exp(-0.1 * (i - 20) * (i - 20)));  // 峰值在i=20
        signal[i] += (uint32_t)(3 * exp(-0.1 * (i - 50) * (i - 50)));  // 峰值在i=50
        signal[i] += (uint32_t)(7 * exp(-0.1 * (i - 80) * (i - 80)));  // 峰值在i=80
        
        // 加入平台峰值
        if (i >= 30 && i <= 35) {
            signal[i] = 4;  // 平台峰值
        }
        
        // 加入隨機噪聲
        signal[i] += (rand() % 2);
    }
    
    return signal;
}

// 將信號保存到CSV文件
void save_signal_to_csv(const uint32_t* signal, int size, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        return;
    }
    
    for (int i = 0; i < size; i++) {
        fprintf(file, "%d,%u\n", i, signal[i]);
    }
    
    fclose(file);
}

// 將峰值結果保存到CSV文件
void save_peaks_to_csv(const uint32_t* signal, const PeakResults* results, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
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
    // 生成測試信號
    int signal_size;
    uint32_t* signal = generate_test_signal(&signal_size);
    save_signal_to_csv(signal, signal_size, "signal_mtk.csv");
    
    // 1. 測試distance參數
    printf("Testing distance parameter:\n");
    PeakResults* results1 = find_peaks(signal, signal_size, 
                                      15, 0, 0, 0, 0, 0, 0, 0, 50);
    
    printf("Found %d peaks using distance=15\n", results1->size);
    for (int i = 0; i < results1->size; i++) {
        printf("  Peak #%d: index=%d, height=%u\n", 
               i+1, results1->indices[i], signal[results1->indices[i]]);
    }
    save_peaks_to_csv(signal, results1, "peaks_distance_mtk.csv");
    free_peak_results(results1);
    
    // 2. 測試prominence參數
    printf("\nTesting prominence parameter:\n");
    PeakResults* results2 = find_peaks(signal, signal_size, 
                                      0, 3, 0, 0, 0, 0, 0, 0, 50);
    
    printf("Found %d peaks using min_prominence=3\n", results2->size);
    for (int i = 0; i < results2->size; i++) {
        printf("  Peak #%d: index=%d, height=%u, prominence=%u\n", 
               i+1, results2->indices[i], signal[results2->indices[i]], 
               results2->prominences[i]);
    }
    save_peaks_to_csv(signal, results2, "peaks_prominence_mtk.csv");
    free_peak_results(results2);
    
    // 3. 測試width參數
    printf("\nTesting width parameter:\n");
    PeakResults* results3 = find_peaks(signal, signal_size, 
                                      0, 0, 0, 5, 0, 0, 0, 0, 50);
    
    printf("Found %d peaks using min_width=5\n", results3->size);
    for (int i = 0; i < results3->size; i++) {
        printf("  Peak #%d: index=%d, height=%u, width=%u\n", 
               i+1, results3->indices[i], signal[results3->indices[i]], 
               results3->widths[i]);
    }
    save_peaks_to_csv(signal, results3, "peaks_width_mtk.csv");
    free_peak_results(results3);
    
    // 4. 測試plateau_size參數
    printf("\nTesting plateau size parameter:\n");
    PeakResults* results4 = find_peaks(signal, signal_size, 
                                      0, 0, 0, 0, 0, 3, 0, 0, 50);
    
    printf("Found %d peaks using min_plateau_size=3\n", results4->size);
    for (int i = 0; i < results4->size; i++) {
        printf("  Peak #%d: index=%d, height=%u, plateau_size=%d\n", 
               i+1, results4->indices[i], signal[results4->indices[i]], 
               results4->plateau_sizes[i]);
    }
    save_peaks_to_csv(signal, results4, "peaks_plateau_mtk.csv");
    free_peak_results(results4);
    
    // 5. 測試組合參數
    printf("\nTesting combined parameters:\n");
    PeakResults* results5 = find_peaks(signal, signal_size, 
                                      10, 2, 0, 3, 0, 0, 0, 0, 50);
    
    printf("Found %d peaks using combined parameters\n", results5->size);
    for (int i = 0; i < results5->size; i++) {
        printf("  Peak #%d: index=%d, height=%u, prominence=%u, width=%u\n", 
               i+1, results5->indices[i], signal[results5->indices[i]], 
               results5->prominences[i], results5->widths[i]);
    }
    save_peaks_to_csv(signal, results5, "peaks_combined_mtk.csv");
    free_peak_results(results5);
    
    // 釋放信號記憶體
    free_ctrl_buffer(signal);
    
    return 0;
} 