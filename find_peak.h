/**
 * find_peak.h - C版本的峰值檢測函數
 * 轉換自Python的scipy.signal模塊中的find_peaks功能
 */

#ifndef FIND_PEAK_H
#define FIND_PEAK_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <string.h>

// 定義返回的峰值屬性結構
typedef struct {
    double* peak_heights;       // 峰值高度
    double* left_thresholds;    // 左側閾值
    double* right_thresholds;   // 右側閾值
    double* prominences;        // 顯著性
    size_t* left_bases;         // 左側基底
    size_t* right_bases;        // 右側基底
    double* widths;             // 寬度
    double* width_heights;      // 寬度高度
    double* left_ips;           // 左側交叉點
    double* right_ips;          // 右側交叉點
    size_t* plateau_sizes;      // 平台大小
    size_t* left_edges;         // 左側邊緣
    size_t* right_edges;        // 右側邊緣
    size_t count;               // 峰值數量
} PeakProperties;

// 主要的峰值查找函數
/**
 * 查找信號中的峰值
 * @param x 輸入信號數組
 * @param n 信號長度
 * @param height_min 最小高度 (可為NULL)
 * @param height_max 最大高度 (可為NULL)
 * @param threshold_min 最小閾值 (可為NULL)
 * @param threshold_max 最大閾值 (可為NULL)
 * @param distance 峰值間最小距離 (可為0表示不考慮)
 * @param prominence_min 最小顯著性 (可為NULL)
 * @param prominence_max 最大顯著性 (可為NULL)
 * @param width_min 最小寬度 (可為NULL)
 * @param width_max 最大寬度 (可為NULL)
 * @param wlen 用於計算顯著性的窗口大小 (可為0表示使用全部信號)
 * @param rel_height 用於計算寬度的相對高度 (0-1之間)
 * @param plateau_size_min 最小平台大小 (可為NULL)
 * @param plateau_size_max 最大平台大小 (可為NULL)
 * @param peaks 輸出的峰值索引數組 (需要使用free釋放)
 * @param properties 輸出的峰值屬性 (需要使用free_peak_properties釋放)
 * @return 峰值數量
 */
size_t find_peaks(const double* x, size_t n,
                 const double* height_min, const double* height_max,
                 const double* threshold_min, const double* threshold_max,
                 size_t distance,
                 const double* prominence_min, const double* prominence_max,
                 const double* width_min, const double* width_max,
                 size_t wlen, double rel_height,
                 const size_t* plateau_size_min, const size_t* plateau_size_max,
                 size_t** peaks, PeakProperties* properties);

/**
 * 使用小波變換查找信號中的峰值
 * @param vector 輸入信號數組
 * @param n 信號長度
 * @param widths 小波寬度數組
 * @param n_widths 小波寬度數量
 * @param min_snr 最小信噪比
 * @param noise_perc 噪聲百分比
 * @param peaks 輸出的峰值索引數組 (需要使用free釋放)
 * @return 峰值數量
 */
size_t find_peaks_cwt(const double* vector, size_t n, 
                     const double* widths, size_t n_widths,
                     double min_snr, double noise_perc,
                     size_t** peaks);

/**
 * 釋放峰值屬性結構中的內存
 * @param props 要釋放的峰值屬性結構
 */
void free_peak_properties(PeakProperties* props);

/**
 * 計算峰值的顯著性
 * @param x 輸入信號數組
 * @param n 信號長度
 * @param peaks 峰值索引數組
 * @param n_peaks 峰值數量
 * @param wlen 窗口長度 (可為0表示使用全部信號)
 * @param prominences 輸出的顯著性數組 (需要預先分配)
 * @param left_bases 輸出的左側基底索引 (需要預先分配)
 * @param right_bases 輸出的右側基底索引 (需要預先分配)
 */
void peak_prominences(const double* x, size_t n,
                     const size_t* peaks, size_t n_peaks,
                     size_t wlen,
                     double* prominences, 
                     size_t* left_bases, 
                     size_t* right_bases);

/**
 * 計算峰值的寬度
 * @param x 輸入信號數組
 * @param n 信號長度
 * @param peaks 峰值索引數組
 * @param n_peaks 峰值數量
 * @param rel_height 相對高度 (0-1之間)
 * @param prominences 峰值顯著性數組
 * @param left_bases 左側基底索引數組
 * @param right_bases 右側基底索引數組
 * @param widths 輸出的寬度數組 (需要預先分配)
 * @param width_heights 輸出的寬度高度數組 (需要預先分配)
 * @param left_ips 輸出的左側交叉點數組 (需要預先分配)
 * @param right_ips 輸出的右側交叉點數組 (需要預先分配)
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