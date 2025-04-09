/**
 * find_peak.c - C版本的峰值檢測函數的實現
 * 轉換自Python的scipy.signal模塊中的find_peaks功能
 */

#include "find_peak.h"

// 初始化峰值屬性結構
static void init_peak_properties(PeakProperties* props, size_t count) {
    props->count = count;
    props->peak_heights = NULL;
    props->left_thresholds = NULL;
    props->right_thresholds = NULL;
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

// 釋放峰值屬性結構中的內存
void free_peak_properties(PeakProperties* props) {
    if (props == NULL) return;
    
    free(props->peak_heights);
    free(props->left_thresholds);
    free(props->right_thresholds);
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
    
    // 重置所有指針
    init_peak_properties(props, 0);
}

// 排序數組並返回排序後的索引
static void argsort(const double* array, size_t n, size_t* indices) {
    // 初始化索引數組
    for (size_t i = 0; i < n; i++) {
        indices[i] = i;
    }
    
    // 使用插入排序對索引進行排序
    for (size_t i = 1; i < n; i++) {
        size_t j = i;
        while (j > 0 && array[indices[j-1]] > array[indices[j]]) {
            // 交換索引
            size_t temp = indices[j];
            indices[j] = indices[j-1];
            indices[j-1] = temp;
            j--;
        }
    }
}

// 找出局部極大值
static size_t _local_maxima_1d(const double* x, size_t n, 
                              size_t** peaks, 
                              size_t** left_edges, 
                              size_t** right_edges) {
    // 分配臨時數組來標記局部極大值
    bool* is_max = (bool*)calloc(n, sizeof(bool));
    if (!is_max) {
        return 0; // 內存分配失敗
    }
    
    // 找出所有局部極大值
    size_t max_count = 0;
    
    // 特殊處理第一個點
    if (n > 1 && x[0] > x[1]) {
        is_max[0] = true;
        max_count++;
    }
    
    // 處理中間的點
    for (size_t i = 1; i < n - 1; i++) {
        if (x[i] > x[i-1] && x[i] >= x[i+1]) {
            is_max[i] = true;
            max_count++;
        }
    }
    
    // 特殊處理最後一個點
    if (n > 1 && x[n-1] > x[n-2]) {
        is_max[n-1] = true;
        max_count++;
    }
    
    // 分配結果數組
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
        return 0; // 內存分配失敗
    }
    
    // 填充結果數組
    size_t count = 0;
    
    // 遍歷找出平台邊緣
    for (size_t i = 0; i < n; i++) {
        if (is_max[i]) {
            // 找到左邊緣
            size_t left = i;
            while (left > 0 && x[left-1] == x[i]) {
                left--;
            }
            
            // 找到右邊緣
            size_t right = i;
            while (right < n - 1 && x[right+1] == x[i]) {
                right++;
            }
            
            // 計算峰值位置(對於平台取中間位置)
            size_t peak = left + (right - left) / 2;
            
            (*peaks)[count] = peak;
            (*left_edges)[count] = left;
            (*right_edges)[count] = right;
            count++;
            
            // 跳過已經處理過的平台
            i = right;
        }
    }
    
    free(is_max);
    return count;
}

// 計算數據的百分位數
static double scoreatpercentile(const double* data, size_t n, double per) {
    if (n == 0) return 0.0;
    if (n == 1) return data[0];
    
    // 複製數據以便排序
    double* sorted = (double*)malloc(n * sizeof(double));
    if (!sorted) return 0.0;
    
    memcpy(sorted, data, n * sizeof(double));
    
    // 簡單的插入排序
    for (size_t i = 1; i < n; i++) {
        double key = sorted[i];
        int j = i - 1;
        while (j >= 0 && sorted[j] > key) {
            sorted[j + 1] = sorted[j];
            j--;
        }
        sorted[j + 1] = key;
    }
    
    // 計算索引
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

// 根據屬性選擇峰值
static bool* _select_by_property(const double* property, size_t n,
                               const double* pmin, const double* pmax) {
    bool* keep = (bool*)malloc(n * sizeof(bool));
    if (!keep) return NULL;
    
    // 初始化全部保留
    for (size_t i = 0; i < n; i++) {
        keep[i] = true;
    }
    
    // 應用最小閾值
    if (pmin != NULL) {
        for (size_t i = 0; i < n; i++) {
            if (property[i] < *pmin) {
                keep[i] = false;
            }
        }
    }
    
    // 應用最大閾值
    if (pmax != NULL) {
        for (size_t i = 0; i < n; i++) {
            if (property[i] > *pmax) {
                keep[i] = false;
            }
        }
    }
    
    return keep;
}

// 根據峰值閾值選擇峰值
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
    
    // 初始化全部保留
    for (size_t i = 0; i < n; i++) {
        keep[i] = true;
    }
    
    // 計算每個峰值的左右閾值
    for (size_t i = 0; i < n; i++) {
        size_t idx = peaks[i];
        double left_threshold = 0, right_threshold = 0;
        
        // 確保索引在有效範圍內
        if (idx > 0) {
            left_threshold = x[idx] - x[idx-1];
        }
        
        if (idx < n - 1) {
            right_threshold = x[idx] - x[idx+1];
        }
        
        (*left_thresholds)[i] = left_threshold;
        (*right_thresholds)[i] = right_threshold;
        
        // 計算最小閾值
        double min_threshold = (left_threshold < right_threshold) ? left_threshold : right_threshold;
        
        // 計算最大閾值
        double max_threshold = (left_threshold > right_threshold) ? left_threshold : right_threshold;
        
        // 應用最小閾值條件
        if (tmin != NULL && min_threshold < *tmin) {
            keep[i] = false;
        }
        
        // 應用最大閾值條件
        if (tmax != NULL && max_threshold > *tmax) {
            keep[i] = false;
        }
    }
    
    return keep;
}

// 根據峰值間距選擇峰值
static bool* _select_by_peak_distance(const size_t* peaks, size_t n, const double* values, size_t distance) {
    if (distance <= 1 || n == 0) {
        // 如果沒有距離限制或沒有峰值，保留所有峰值
        bool* keep = (bool*)malloc(n * sizeof(bool));
        if (!keep) return NULL;
        
        for (size_t i = 0; i < n; i++) {
            keep[i] = true;
        }
        return keep;
    }
    
    // 首先按照峰值高度排序
    size_t* sorted_indices = (size_t*)malloc(n * sizeof(size_t));
    if (!sorted_indices) return NULL;
    
    argsort(values, n, sorted_indices);
    
    // 創建結果數組，初始全部保留
    bool* keep = (bool*)malloc(n * sizeof(bool));
    if (!keep) {
        free(sorted_indices);
        return NULL;
    }
    
    for (size_t i = 0; i < n; i++) {
        keep[i] = true;
    }
    
    // 按照峰值高度從高到低遍歷
    for (size_t i = n; i > 0; i--) {
        size_t idx = sorted_indices[i-1];
        
        if (keep[idx]) {
            // 標記距離太近的峰值為不保留
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

// 計算峰值的顯著性
void peak_prominences(const double* x, size_t n,
                     const size_t* peaks, size_t n_peaks,
                     size_t wlen,
                     double* prominences, 
                     size_t* left_bases, 
                     size_t* right_bases) {
    // 如果沒有峰值，直接返回
    if (n_peaks == 0) return;
    
    // 如果未指定窗口大小，使用整個信號
    size_t effective_wlen;
    if (wlen == 0 || wlen > n) {
        effective_wlen = n;
    } else {
        // 確保窗口大小為奇數
        effective_wlen = wlen;
        if (effective_wlen % 2 == 0) {
            effective_wlen += 1;
        }
    }
    
    // 處理每個峰值
    for (size_t i = 0; i < n_peaks; i++) {
        size_t peak_idx = peaks[i];
        double peak_height = x[peak_idx];
        
        // 定義窗口範圍
        size_t window_start, window_end;
        
        if (effective_wlen < n) {
            // 使用指定窗口大小
            size_t half_wlen = effective_wlen / 2;
            window_start = (peak_idx > half_wlen) ? peak_idx - half_wlen : 0;
            window_end = (peak_idx + half_wlen < n) ? peak_idx + half_wlen : n - 1;
        } else {
            // 使用整個信號
            window_start = 0;
            window_end = n - 1;
        }
        
        // 在窗口左側找最小值
        size_t left_min_idx = peak_idx;
        double left_min = peak_height;
        bool found_left_crossing = false;
        
        for (size_t j = peak_idx; j > window_start; j--) {
            // 如果找到更高的峰值，停止搜索
            if (x[j] > peak_height) {
                found_left_crossing = true;
                break;
            }
            
            // 更新最小值
            if (x[j] < left_min) {
                left_min = x[j];
                left_min_idx = j;
            }
        }
        
        // 如果窗口太小，檢查窗口起點
        if (!found_left_crossing && window_start > 0) {
            if (x[window_start] < left_min) {
                left_min = x[window_start];
                left_min_idx = window_start;
            }
        }
        
        // 在窗口右側找最小值
        size_t right_min_idx = peak_idx;
        double right_min = peak_height;
        bool found_right_crossing = false;
        
        for (size_t j = peak_idx; j < window_end; j++) {
            // 如果找到更高的峰值，停止搜索
            if (x[j] > peak_height) {
                found_right_crossing = true;
                break;
            }
            
            // 更新最小值
            if (x[j] < right_min) {
                right_min = x[j];
                right_min_idx = j;
            }
        }
        
        // 如果窗口太小，檢查窗口終點
        if (!found_right_crossing && window_end < n - 1) {
            if (x[window_end] < right_min) {
                right_min = x[window_end];
                right_min_idx = window_end;
            }
        }
        
        // 確定峰值的基底高度和顯著性
        double base_height = (left_min > right_min) ? left_min : right_min;
        prominences[i] = peak_height - base_height;
        
        // 保存基底位置
        left_bases[i] = left_min_idx;
        right_bases[i] = right_min_idx;
    }
}

// 線性插值找出交叉點
static double _interpolate_intersection(double x1, double y1, double x2, double y2, double target_y) {
    // 如果兩點的y值相同，無法插值，返回其中一點
    if (y1 == y2) {
        return x1;
    }
    
    // 線性插值計算交叉點
    return x1 + (x2 - x1) * (target_y - y1) / (y2 - y1);
}

// 計算峰值的寬度
void peak_widths(const double* x, size_t n,
                const size_t* peaks, size_t n_peaks,
                double rel_height,
                const double* prominences,
                const size_t* left_bases,
                const size_t* right_bases,
                double* widths,
                double* width_heights,
                double* left_ips,
                double* right_ips) {
    // 如果沒有峰值，直接返回
    if (n_peaks == 0) return;
    
    // 計算每個峰值的寬度
    for (size_t i = 0; i < n_peaks; i++) {
        size_t peak_idx = peaks[i];
        double prominence = prominences[i];
        double peak_height = x[peak_idx];
        
        // 計算評估高度
        double evaluation_height = peak_height - prominence * rel_height;
        width_heights[i] = evaluation_height;
        
        // 找到交叉點位置
        size_t j_left = peak_idx;
        size_t j_right = peak_idx;
        
        // 向左尋找交叉點
        while (j_left > 0 && x[j_left] > evaluation_height) {
            j_left--;
        }
        
        // 向右尋找交叉點
        while (j_right < n - 1 && x[j_right] > evaluation_height) {
            j_right++;
        }
        
        // 如果達到了信號邊界，使用邊界值
        if (j_left == 0 && x[j_left] > evaluation_height) {
            left_ips[i] = 0.0;
        } else {
            // 通過線性插值計算精確的交叉點位置
            size_t idx_left = j_left;
            size_t idx_right = j_left + 1;
            
            // 確保不越界
            if (idx_right >= n) {
                idx_right = n - 1;
            }
            
            left_ips[i] = _interpolate_intersection(idx_left, x[idx_left], idx_right, x[idx_right], evaluation_height);
        }
        
        // 如果達到了信號邊界，使用邊界值
        if (j_right == n - 1 && x[j_right] > evaluation_height) {
            right_ips[i] = n - 1.0;
        } else {
            // 通過線性插值計算精確的交叉點位置
            size_t idx_left = j_right - 1;
            size_t idx_right = j_right;
            
            // 確保不越界
            if (idx_left >= n) {
                idx_left = n - 1;
            }
            
            right_ips[i] = _interpolate_intersection(idx_left, x[idx_left], idx_right, x[idx_right], evaluation_height);
        }
        
        // 計算寬度
        widths[i] = right_ips[i] - left_ips[i];
    }
}

// 主要的峰值查找函數
size_t find_peaks(const double* x, size_t n,
                 const double* height_min, const double* height_max,
                 const double* threshold_min, const double* threshold_max,
                 size_t distance,
                 const double* prominence_min, const double* prominence_max,
                 const double* width_min, const double* width_max,
                 size_t wlen, double rel_height,
                 const size_t* plateau_size_min, const size_t* plateau_size_max,
                 size_t** peaks, PeakProperties* properties) {
    // 初始化輸出
    *peaks = NULL;
    init_peak_properties(properties, 0);
    
    if (n == 0) {
        return 0;
    }
    
    // 檢查距離參數
    if (distance > 0 && distance < 1) {
        fprintf(stderr, "Error: distance must be greater or equal to 1\n");
        return 0;
    }
    
    // 找出所有的局部極大值
    size_t* all_peaks = NULL;
    size_t* left_edges = NULL;
    size_t* right_edges = NULL;
    size_t peak_count = _local_maxima_1d(x, n, &all_peaks, &left_edges, &right_edges);
    
    if (peak_count == 0) {
        return 0; // 沒有找到峰值
    }
    
    // 保存所有峰值的副本，用於後續篩選
    size_t* current_peaks = (size_t*)malloc(peak_count * sizeof(size_t));
    if (!current_peaks) {
        free(all_peaks);
        free(left_edges);
        free(right_edges);
        return 0;
    }
    
    memcpy(current_peaks, all_peaks, peak_count * sizeof(size_t));
    size_t current_count = peak_count;
    
    // 評估平台大小
    if (plateau_size_min != NULL || plateau_size_max != NULL) {
        // 計算平台大小
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
        
        // 篩選符合平台大小條件的峰值
        bool* keep = _select_by_property(plateau_sizes, current_count, 
                                        plateau_size_min, plateau_size_max);
        if (!keep) {
            free(all_peaks);
            free(left_edges);
            free(right_edges);
            free(current_peaks);
            free(plateau_sizes);
            return 0;
        }
        
        // 更新峰值列表
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
        
        // 保存平台相關的屬性
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
    
    // 評估高度條件
    if (height_min != NULL || height_max != NULL) {
        // 獲取峰值高度
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
        
        // 篩選符合高度條件的峰值
        bool* keep = _select_by_property(peak_heights, current_count, height_min, height_max);
        if (!keep) {
            free(all_peaks);
            free(left_edges);
            free(right_edges);
            free(current_peaks);
            free(peak_heights);
            free_peak_properties(properties);
            return 0;
        }
        
        // 更新峰值列表
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
        
        // 保存高度屬性
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
    
    // 評估閾值條件
    if (threshold_min != NULL || threshold_max != NULL) {
        double* left_thresholds = NULL;
        double* right_thresholds = NULL;
        
        // 篩選符合閾值條件的峰值
        bool* keep = _select_by_peak_threshold(x, current_peaks, current_count, 
                                             threshold_min, threshold_max,
                                             &left_thresholds, &right_thresholds);
        if (!keep) {
            free(all_peaks);
            free(left_edges);
            free(right_edges);
            free(current_peaks);
            free_peak_properties(properties);
            return 0;
        }
        
        // 更新峰值列表
        size_t new_count = 0;
        for (size_t i = 0; i < current_count; i++) {
            if (keep[i]) {
                current_peaks[new_count] = current_peaks[i];
                left_thresholds[new_count] = left_thresholds[i];
                right_thresholds[new_count] = right_thresholds[i];
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
        
        // 保存閾值屬性
        if (current_count > 0) {
            properties->left_thresholds = (double*)malloc(current_count * sizeof(double));
            properties->right_thresholds = (double*)malloc(current_count * sizeof(double));
            
            if (!properties->left_thresholds || !properties->right_thresholds) {
                free(all_peaks);
                free(left_edges);
                free(right_edges);
                free(current_peaks);
                free(left_thresholds);
                free(right_thresholds);
                free(keep);
                free_peak_properties(properties);
                return 0;
            }
            
            memcpy(properties->left_thresholds, left_thresholds, current_count * sizeof(double));
            memcpy(properties->right_thresholds, right_thresholds, current_count * sizeof(double));
        }
        
        free(keep);
        free(left_thresholds);
        free(right_thresholds);
    }
    
    // 評估距離條件
    if (distance > 1) {
        // 獲取峰值高度(如果尚未計算)
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
        
        // 篩選符合距離條件的峰值
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
        
        // 更新峰值列表
        size_t new_count = 0;
        for (size_t i = 0; i < current_count; i++) {
            if (keep[i]) {
                current_peaks[new_count] = current_peaks[i];
                if (properties->peak_heights) {
                    properties->peak_heights[new_count] = properties->peak_heights[i];
                }
                if (properties->left_thresholds) {
                    properties->left_thresholds[new_count] = properties->left_thresholds[i];
                    properties->right_thresholds[new_count] = properties->right_thresholds[i];
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
        
        // 如果新生成了peak_heights，則釋放它
        if (!properties->peak_heights) {
            free(peak_heights);
        }
        
        free(keep);
    }
    
    // 評估顯著性條件
    if (prominence_min != NULL || prominence_max != NULL || width_min != NULL || width_max != NULL) {
        // 計算顯著性
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
        
        // 篩選符合顯著性條件的峰值
        if (prominence_min != NULL || prominence_max != NULL) {
            bool* keep = _select_by_property(prominences, current_count, prominence_min, prominence_max);
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
            
            // 更新峰值列表
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
                    if (properties->left_thresholds) {
                        properties->left_thresholds[new_count] = properties->left_thresholds[i];
                        properties->right_thresholds[new_count] = properties->right_thresholds[i];
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
        
        // 保存顯著性屬性
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
        
        // 評估寬度條件
        if (width_min != NULL || width_max != NULL) {
            // 計算寬度
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
            
            // 篩選符合寬度條件的峰值
            bool* keep = _select_by_property(widths, current_count, width_min, width_max);
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
            
            // 更新峰值列表
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
                    if (properties->left_thresholds) {
                        properties->left_thresholds[new_count] = properties->left_thresholds[i];
                        properties->right_thresholds[new_count] = properties->right_thresholds[i];
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
            
            // 保存寬度屬性
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
    
    // 返回最終的峰值列表
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
    
    // 釋放臨時數組
    free(all_peaks);
    free(left_edges);
    free(right_edges);
    free(current_peaks);
    
    return current_count;
}

// Ricker小波(墨西哥帽小波)
static void _ricker(double* wavelet, size_t n, double width) {
    double a = (n - 1) / 2.0;
    double wsq = width * width;
    
    for (size_t i = 0; i < n; i++) {
        double x = i - a;
        double xsq = x * x;
        double mod = 1.0 - xsq / wsq;
        double gauss = exp(-xsq / (2.0 * wsq));
        wavelet[i] = mod * gauss;
    }
}

// 小波轉換函數
static void _cwt(const double* data, size_t n, double** out, 
                const double* widths, size_t n_widths) {
    // 分配輸出數組
    for (size_t i = 0; i < n_widths; i++) {
        out[i] = (double*)calloc(n, sizeof(double));
        if (!out[i]) {
            return;
        }
    }
    
    // 對每個寬度執行小波轉換
    for (size_t i = 0; i < n_widths; i++) {
        double width = widths[i];
        
        // 計算小波大小
        size_t wavelet_size = (size_t)(10 * width);
        if (wavelet_size % 2 == 0) {
            wavelet_size += 1;
        }
        
        // 生成小波
        double* wavelet = (double*)malloc(wavelet_size * sizeof(double));
        if (!wavelet) {
            return;
        }
        
        _ricker(wavelet, wavelet_size, width);
        
        // 標準化小波
        double sum = 0.0;
        for (size_t j = 0; j < wavelet_size; j++) {
            sum += fabs(wavelet[j]);
        }
        
        for (size_t j = 0; j < wavelet_size; j++) {
            wavelet[j] /= sum;
        }
        
        // 卷積
        int half_wavelet_size = (int)(wavelet_size / 2);
        for (size_t j = 0; j < n; j++) {
            for (int k = -half_wavelet_size; k <= half_wavelet_size; k++) {
                int idx = (int)j + k;
                if (idx >= 0 && idx < (int)n) {
                    out[i][j] += data[idx] * wavelet[k + half_wavelet_size];
                }
            }
        }
        
        free(wavelet);
    }
}

// 使用小波變換查找峰值
size_t find_peaks_cwt(const double* vector, size_t n, 
                     const double* widths, size_t n_widths,
                     double min_snr, double noise_perc,
                     size_t** peaks) {
    *peaks = NULL;
    
    if (n == 0 || n_widths == 0) {
        return 0;
    }
    
    // 分配內存用於CWT結果
    double** cwt_matrix = (double**)malloc(n_widths * sizeof(double*));
    if (!cwt_matrix) {
        return 0;
    }
    
    // 執行小波變換
    _cwt(vector, n, cwt_matrix, widths, n_widths);
    
    // 檢查內存分配是否成功
    for (size_t i = 0; i < n_widths; i++) {
        if (!cwt_matrix[i]) {
            for (size_t j = 0; j < i; j++) {
                free(cwt_matrix[j]);
            }
            free(cwt_matrix);
            return 0;
        }
    }
    
    // 使用最小尺度下的訊號作為基準
    const double* row_one = cwt_matrix[0];
    double* noises = (double*)malloc(n * sizeof(double));
    if (!noises) {
        for (size_t i = 0; i < n_widths; i++) {
            free(cwt_matrix[i]);
        }
        free(cwt_matrix);
        return 0;
    }
    
    // 計算噪聲閾值
    size_t window_size = n / 20;
    if (window_size < 1) window_size = 1;
    
    size_t half_window = window_size / 2;
    
    for (size_t i = 0; i < n; i++) {
        size_t window_start = i > half_window ? i - half_window : 0;
        size_t window_end = i + half_window < n ? i + half_window : n;
        size_t window_width = window_end - window_start;
        
        // 提取窗口數據
        double* window_data = (double*)malloc(window_width * sizeof(double));
        if (!window_data) {
            for (size_t j = 0; j < n_widths; j++) {
                free(cwt_matrix[j]);
            }
            free(cwt_matrix);
            free(noises);
            return 0;
        }
        
        for (size_t j = 0; j < window_width; j++) {
            window_data[j] = row_one[window_start + j];
        }
        
        noises[i] = scoreatpercentile(window_data, window_width, noise_perc);
        free(window_data);
    }
    
    // 找出所有可能的峰值(在最小尺度下)
    size_t* max_locs = NULL;
    size_t* left_edges = NULL;
    size_t* right_edges = NULL;
    size_t peak_count = _local_maxima_1d(row_one, n, &max_locs, &left_edges, &right_edges);
    
    // 如果沒有找到峰值，釋放內存並返回
    if (peak_count == 0) {
        for (size_t i = 0; i < n_widths; i++) {
            free(cwt_matrix[i]);
        }
        free(cwt_matrix);
        free(noises);
        return 0;
    }
    
    // 篩選具有足夠信噪比的峰值
    bool* keep = (bool*)malloc(peak_count * sizeof(bool));
    if (!keep) {
        for (size_t i = 0; i < n_widths; i++) {
            free(cwt_matrix[i]);
        }
        free(cwt_matrix);
        free(noises);
        free(max_locs);
        free(left_edges);
        free(right_edges);
        return 0;
    }
    
    size_t valid_count = 0;
    for (size_t i = 0; i < peak_count; i++) {
        size_t loc = max_locs[i];
        double snr = fabs(row_one[loc] / noises[loc]);
        keep[i] = (snr >= min_snr);
        if (keep[i]) {
            valid_count++;
        }
    }
    
    // 分配最終結果數組
    *peaks = (size_t*)malloc(valid_count * sizeof(size_t));
    if (!*peaks) {
        for (size_t i = 0; i < n_widths; i++) {
            free(cwt_matrix[i]);
        }
        free(cwt_matrix);
        free(noises);
        free(max_locs);
        free(left_edges);
        free(right_edges);
        free(keep);
        return 0;
    }
    
    // 填充結果數組
    size_t result_idx = 0;
    for (size_t i = 0; i < peak_count; i++) {
        if (keep[i]) {
            (*peaks)[result_idx++] = max_locs[i];
        }
    }
    
    // 排序峰值位置
    for (size_t i = 1; i < valid_count; i++) {
        size_t key = (*peaks)[i];
        int j = i - 1;
        while (j >= 0 && (*peaks)[j] > key) {
            (*peaks)[j + 1] = (*peaks)[j];
            j--;
        }
        (*peaks)[j + 1] = key;
    }
    
    // 釋放所有臨時內存
    for (size_t i = 0; i < n_widths; i++) {
        free(cwt_matrix[i]);
    }
    free(cwt_matrix);
    free(noises);
    free(max_locs);
    free(left_edges);
    free(right_edges);
    free(keep);
    
    return valid_count;
} 