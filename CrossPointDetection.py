from typing import List, Tuple
import matplotlib.pyplot as plt
import random
import numpy as np
import os
import pickle
from scipy.stats import pearsonr  # 引入pearsonr用於計算相關係數

def find_intersections(lines: List[List[Tuple[float, float]]]) -> List[Tuple[float, float]]:
    def ccw(A, B, C):
        return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])

    def intersect(A, B, C, D):
        return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)

    # 計算線段的斜率
    def calculate_slope(p1, p2):
        if p2[0] - p1[0] == 0:  # 避免除以零
            return float('inf') if p2[1] > p1[1] else float('-inf')
        return (p2[1] - p1[1]) / (p2[0] - p1[0])

    intersections = []
    intersection_info = []  # 存儲交點及相關信息 (交點, 線段1斜率, 線段2斜率, 線1索引, 線2索引, 線段1索引, 線段2索引)
    
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            for k in range(len(lines[i]) - 1):
                for l in range(len(lines[j]) - 1):
                    A, B = lines[i][k], lines[i][k + 1]
                    C, D = lines[j][l], lines[j][l + 1]
                    if intersect(A, B, C, D):
                        # 計算斜率
                        slope1 = calculate_slope(A, B)
                        slope2 = calculate_slope(C, D)
                        
                        # 計算交點
                        denom = (A[0] - B[0]) * (C[1] - D[1]) - (A[1] - B[1]) * (C[0] - D[0])
                        if denom == 0:
                            continue  # 線段平行
                        x = ((A[0] * B[1] - A[1] * B[0]) * (C[0] - D[0]) - (A[0] - B[0]) * (C[0] * D[1] - C[1] * D[0])) / denom
                        y = ((A[0] * B[1] - A[1] * B[0]) * (C[1] - D[1]) - (A[1] - B[1]) * (C[0] * D[1] - C[1] * D[0])) / denom
                        
                        intersection = (x, y)
                        intersections.append(intersection)
                        intersection_info.append((intersection, slope1, slope2, i, j, k, l))
    
    return intersections, intersection_info

# 計算交點附近更大範圍的斜率和相關係數
def calculate_wider_slope(line, segment_idx, window_size=3):
    """計算線段附近更大範圍的平均斜率和時間相關性"""
    # 確保不超出線段範圍
    start_idx = max(0, segment_idx - window_size)
    end_idx = min(len(line) - 1, segment_idx + window_size + 1)
    
    if end_idx <= start_idx + 1:  # 確保至少有兩個點可計算斜率
        return 0, 0, 0
    
    # 取更大範圍的起點和終點
    start_point = line[start_idx]
    end_point = line[end_idx]
    
    # 計算整體斜率
    if end_point[0] - start_point[0] == 0:  # 避免除以零
        return float('inf') if end_point[1] > start_point[1] else float('-inf'), 0, 0
    
    overall_slope = (end_point[1] - start_point[1]) / (end_point[0] - start_point[0])
    
    # 提取範圍內的x和y值
    x_values = [p[0] for p in line[start_idx:end_idx+1]]
    y_values = [p[1] for p in line[start_idx:end_idx+1]]
    
    # 計算時間和信號強度的相關係數
    # 如果點數太少則可能無法計算相關係數
    if len(x_values) > 2:
        try:
            correlation, _ = pearsonr(x_values, y_values)
        except:
            correlation = 0  # 如果計算失敗，設為零
    else:
        correlation = 0
    
    # 計算中間點的斜率一致性，檢測是否有太多震盪
    slopes = []
    for i in range(start_idx, end_idx):
        if i + 1 <= end_idx:
            p1, p2 = line[i], line[i + 1]
            if p2[0] - p1[0] != 0:
                slopes.append((p2[1] - p1[1]) / (p2[0] - p1[0]))
    
    # 檢查震盪：計算斜率變號的次數
    sign_changes = 0
    for i in range(1, len(slopes)):
        if slopes[i] * slopes[i-1] < 0:  # 斜率變號
            sign_changes += 1
    
    # 返回斜率、震盪信息和相關係數
    return overall_slope, sign_changes, correlation

# 識別「Good Cross Point」，考慮更大範圍的斜率、相關係數和震盪
def identify_good_cross_points(intersection_info, lines, window_size=3, min_slope_diff=0.2, max_oscillations=4, min_correlation_abs=0.3, min_slope_significance=0.1):
    """
    識別「Good Cross Point」，使用多種數學標準進行評估，不依賴特定點的硬編碼規則。
    
    參數:
        intersection_info: 交點的基本信息
        lines: 所有線段數據
        window_size: 計算斜率時向兩側擴展的窗口大小
        min_slope_diff: 斜率差異的最小閾值
        max_oscillations: 允許的最大震盪次數
        min_correlation_abs: 相關係數的最小絕對值
        min_slope_significance: 斜率絕對值的最小閾值，任何小於這個值的斜率會被視為「不顯著」
    
    返回:
        good_cross_points: 良好交叉點列表
        rejected_points: 被拒絕的交叉點列表
    """
    good_cross_points = []
    rejected_points = []
    
    for point, _, _, line1_idx, line2_idx, seg1_idx, seg2_idx in intersection_info:
        # 計算交點附近更大範圍的斜率、震盪和相關係數
        wider_slope1, oscillations1, correlation1 = calculate_wider_slope(lines[line1_idx], seg1_idx, window_size)
        wider_slope2, oscillations2, correlation2 = calculate_wider_slope(lines[line2_idx], seg2_idx, window_size)
        
        # 斜率差異
        slope_diff = abs(wider_slope1 - wider_slope2)
        
        # 檢查斜率是否顯著
        slope1_significant = abs(wider_slope1) >= min_slope_significance
        slope2_significant = abs(wider_slope2) >= min_slope_significance
        both_slopes_significant = slope1_significant and slope2_significant
        
        # 計算關鍵特徵
        slopes_opposite = (wider_slope1 > 0 and wider_slope2 < 0) or (wider_slope1 < 0 and wider_slope2 > 0)
        correlations_opposite = (correlation1 > 0 and correlation2 < 0) or (correlation1 < 0 and correlation2 > 0)
        strong_correlations = abs(correlation1) >= 0.8 and abs(correlation2) >= 0.8
        
        # 基本震盪檢查，但為特殊情況設定彈性標準
        oscillations_acceptable = oscillations1 <= max_oscillations and oscillations2 <= max_oscillations
        
        # 1. 標準交叉點類型：斜率方向相反
        standard_cross = (
            slopes_opposite and 
            slope_diff >= min_slope_diff and 
            both_slopes_significant and
            oscillations_acceptable and
            abs(correlation1) >= min_correlation_abs and 
            abs(correlation2) >= min_correlation_abs
        )
        
        # 2. 相關係數型交叉點：相關係數方向相反且強相關
        correlation_cross = (
            correlations_opposite and
            strong_correlations and
            both_slopes_significant and
            oscillations_acceptable
        )
        
        # 3. 同向強相關交叉點：斜率雖同向但差異顯著，相關係數強且方向一致
        special_cross = (
            not slopes_opposite and
            slope_diff > 0.3 and
            strong_correlations and
            (correlation1 * correlation2 > 0) and  # 相關係數同向
            both_slopes_significant and
            oscillations_acceptable
        )
        
        # 4. 增強型交叉點：一條線顯示強烈趨勢，可對另一條線適當放寬要求
        # 4.1 - 非平衡斜率交叉點：當一條線有顯著斜率且強相關時，可放寬對另一條線的要求
        strong_slope_threshold = 0.7  # 定義強斜率的閾值
        strong_correlation_threshold = 0.9  # 定義強相關的閾值
        
        line1_strong = abs(wider_slope1) >= strong_slope_threshold and abs(correlation1) >= strong_correlation_threshold
        line2_strong = abs(wider_slope2) >= strong_slope_threshold and abs(correlation2) >= strong_correlation_threshold
        
        asymmetric_cross = (
            slopes_opposite and  # 斜率方向相反
            slope_diff >= min_slope_diff and  # 斜率差異顯著
            oscillations_acceptable and  # 震盪在允許範圍內
            (line1_strong or line2_strong)  # 至少一條線有強烈趨勢
        )
        
        # 4.2 - 顯著斜率差異型：當斜率差異非常顯著時，可放寬其他要求
        very_significant_diff = slope_diff >= 0.6  # 定義顯著斜率差異的閾值
        
        slope_diff_cross = (
            very_significant_diff and  # 斜率差異非常顯著
            oscillations_acceptable and  # 震盪在允許範圍內
            (  # 滿足以下條件之一：
                slopes_opposite or  # 斜率方向相反
                (abs(correlation1) >= 0.7 and abs(correlation2) >= 0.7)  # 兩條線都有較強相關性
            )
        )
        
        # 4.3 - 同向斜率但相關性同向且強：同方向的斜率變化且相關係數一致
        same_direction_cross = (
            not slopes_opposite and  # 斜率方向相同
            slope_diff >= 0.5 and  # 斜率差異夠大
            abs(correlation1) >= 0.65 and abs(correlation2) >= 0.65 and  # 兩條線都有強相關
            (correlation1 * correlation2 > 0) and  # 相關係數同向
            oscillations_acceptable  # 震盪在允許範圍內
        )
        
        # 5. 高震盪容忍型交叉點：在其他條件特別優秀時，容忍稍高的震盪
        # 定義震盪容忍程度 - 根據斜率差和相關性強度動態調整
        max_tolerable_oscillations = 6  # 最高可容忍的震盪次數
        
        # 判斷是否具有超強斜率差和極強相關
        exceptional_slope_diff = slope_diff >= 0.9  # 接近或超過1的斜率差
        exceptional_correlation = (abs(correlation1) >= 0.9 or abs(correlation2) >= 0.9)  # 至少一條線有極強相關
        
        high_oscillation_cross = (
            slopes_opposite and  # 斜率方向相反
            exceptional_slope_diff and  # 超強斜率差
            exceptional_correlation and  # 極強相關
            (oscillations1 <= max_tolerable_oscillations and oscillations2 <= max_tolerable_oscillations) and  # 震盪數在容忍範圍內
            (oscillations1 <= max_oscillations or oscillations2 <= max_oscillations)  # 至少一條線的震盪在標準範圍內
        )
        
        # 6. 弱斜率容忍型交叉點：當對面線條極強時，可以容忍較弱的斜率
        ultra_strong_slope = 0.9  # 極強斜率閾值
        ultra_strong_correlation = 0.95  # 極強相關閾值
        
        line1_ultra_strong = abs(wider_slope1) >= ultra_strong_slope and abs(correlation1) >= ultra_strong_correlation
        line2_ultra_strong = abs(wider_slope2) >= ultra_strong_slope and abs(correlation2) >= ultra_strong_correlation
        
        weak_slope_cross = (
            slopes_opposite and  # 斜率方向相反
            (line1_ultra_strong or line2_ultra_strong) and  # 至少一條線呈現極強趨勢
            oscillations1 <= max_tolerable_oscillations and oscillations2 <= max_tolerable_oscillations and  # 震盪在可容忍範圍
            abs(correlation1) >= min_correlation_abs and abs(correlation2) >= min_correlation_abs  # 基本相關要求
        )
        
        # 7. 相關係數優先型交叉點：基於相關係數方向而非斜率方向判斷交叉
        correlation_priority_cross = (
            correlations_opposite and  # 相關係數方向相反
            abs(correlation1) >= 0.8 and abs(correlation2) >= 0.45 and  # 兩條線都有足夠的相關性
            (slope1_significant or slope2_significant) and  # 至少一條線斜率顯著
            slope_diff >= 0.5 and  # 斜率差異足夠大
            oscillations1 <= max_oscillations and oscillations2 <= max_oscillations  # 震盪在基本範圍內
        )
        
        # 8. NEW: 低震盪強相關交叉點：一條線震盪極低時，允許另一條線震盪稍高
        low_oscillation_threshold = 2  # 定義「低震盪」的閾值
        line1_low_oscillation = oscillations1 <= low_oscillation_threshold
        line2_low_oscillation = oscillations2 <= low_oscillation_threshold
        max_allowed_oscillation = 6  # 當另一條線震盪低時，可容忍的最大震盪值
        
        balanced_oscillation_cross = (
            slopes_opposite and  # 斜率方向相反
            both_slopes_significant and  # 兩條線斜率都顯著
            slope_diff >= 0.4 and  # 斜率差異足夠大
            (  # 至少一條線震盪低，且另一條線震盪不超過最大允許值
                (line1_low_oscillation and oscillations2 <= max_allowed_oscillation) or
                (line2_low_oscillation and oscillations1 <= max_allowed_oscillation)
            ) and
            abs(correlation1) >= 0.6 and abs(correlation2) >= 0.6 and  # 兩條線相關性都較強
            correlations_opposite  # 相關係數方向相反
        )
        
        # 9. NEW: 斜率方向相反但一線震盪稍高的交叉點
        opposite_slope_priority = (
            slopes_opposite and  # 斜率方向相反
            both_slopes_significant and  # 兩條線斜率都顯著
            slope_diff >= 0.4 and  # 斜率差異足夠大
            correlations_opposite and  # 相關係數方向相反
            abs(correlation1) >= 0.65 and abs(correlation2) >= 0.65 and  # 較強相關
            (  # 其中一條線震盪在正常範圍內，另一條只稍微超出
                (oscillations1 <= max_oscillations and oscillations2 <= max_oscillations + 2) or
                (oscillations2 <= max_oscillations and oscillations1 <= max_oscillations + 2)
            )
        )
        
        # 決定此點是否為良好交點
        is_good_point = (
            standard_cross or 
            correlation_cross or 
            special_cross or 
            asymmetric_cross or 
            slope_diff_cross or 
            same_direction_cross or
            high_oscillation_cross or
            weak_slope_cross or
            correlation_priority_cross or
            balanced_oscillation_cross or
            opposite_slope_priority
        )
        
        if is_good_point:
            # 為良好交點添加分類原因
            if standard_cross:
                reason = "標準交叉點：斜率方向相反且差異顯著"
            elif correlation_cross:
                reason = "相關係數型交叉點：相關係數方向相反且強相關"
            elif special_cross:
                reason = "同向強相關交叉點：斜率同向但差異顯著，相關係數強且一致"
            elif asymmetric_cross:
                reason = "非平衡交叉點：一條線有強烈趨勢，另一條線有相反方向的變化"
            elif slope_diff_cross:
                reason = "顯著斜率差異型交叉點：斜率差異極為顯著"
            elif same_direction_cross:
                reason = "同步趨勢型交叉點：斜率方向相同但差異大，相關係數同向且強"
            elif high_oscillation_cross:
                reason = "高震盪容忍型交叉點：雖震盪稍高但斜率差異和相關性極佳"
            elif weak_slope_cross:
                reason = "弱斜率容忍型交叉點：雖一線斜率較弱但對面線條呈極強趨勢"
            elif correlation_priority_cross:
                reason = "相關係數優先型交叉點：基於相關係數方向判斷交叉而非斜率"
            elif balanced_oscillation_cross:
                reason = "低震盪強相關交叉點：一條線震盪極低時容忍另一條線震盪稍高"
            elif opposite_slope_priority:
                reason = "斜率優先型交叉點：斜率方向相反且差異明顯，容忍一線震盪稍高"
            else:
                reason = "綜合條件良好的交叉點"
                
            good_cross_points.append((point, wider_slope1, wider_slope2, line1_idx, line2_idx, 
                                    oscillations1, oscillations2, correlation1, correlation2, reason))
        else:
            # 記錄被拒絕的點及原因
            rejection_reason = []
            
            # 根據具體情況給出拒絕原因
            if not slope1_significant:
                rejection_reason.append(f"線1斜率不夠顯著 (|{wider_slope1:.3f}| < {min_slope_significance})")
            if not slope2_significant:
                rejection_reason.append(f"線2斜率不夠顯著 (|{wider_slope2:.3f}| < {min_slope_significance})")
                
            # 檢查斜率方向和差異
            if not slopes_opposite and not (special_cross or same_direction_cross or correlation_priority_cross):
                rejection_reason.append("斜率不是一正一負且不符合特殊交叉條件")
            if slope_diff < min_slope_diff:
                rejection_reason.append(f"斜率差異不足 ({slope_diff:.3f} < {min_slope_diff})")
                
            # 檢查震盪
            if oscillations1 > max_tolerable_oscillations:
                rejection_reason.append(f"線1震盪嚴重過多 ({oscillations1} > {max_tolerable_oscillations})")
            elif oscillations1 > max_oscillations:
                rejection_reason.append(f"線1震盪過多 ({oscillations1} > {max_oscillations})")
                
            if oscillations2 > max_tolerable_oscillations:
                rejection_reason.append(f"線2震盪嚴重過多 ({oscillations2} > {max_tolerable_oscillations})")
            elif oscillations2 > max_oscillations:
                rejection_reason.append(f"線2震盪過多 ({oscillations2} > {max_oscillations})")
                
            # 檢查相關性
            if abs(correlation1) < min_correlation_abs:
                rejection_reason.append(f"線1相關性不足 (|{correlation1:.3f}| < {min_correlation_abs})")
            if abs(correlation2) < min_correlation_abs:
                rejection_reason.append(f"線2相關性不足 (|{correlation2:.3f}| < {min_correlation_abs})")
            
            rejected_points.append((point, wider_slope1, wider_slope2, line1_idx, line2_idx, 
                                oscillations1, oscillations2, correlation1, correlation2, ", ".join(rejection_reason)))
    
    return good_cross_points, rejected_points

# 生成更加多變的時間序列數據
def generate_time_series(num_series: int, num_points: int) -> List[List[Tuple[float, float]]]:
    series_list = []
    time_range = np.linspace(0, 100, num_points)  # 時間從0到100
    
    # 設定隨機種子使結果可重複但與之前不同
    random.seed(987)  # 使用與之前不同的種子
    
    for i in range(num_series):
        # 使用更複雜的信號模型
        # 基本正弦波 + 二次趨勢 + 隨機噪聲
        amplitude = random.uniform(2, 7)  # 更大的振幅範圍
        frequency = random.uniform(0.05, 0.3)  # 不同的頻率範圍
        phase = random.uniform(0, 2*np.pi)  # 隨機相位
        offset = random.uniform(-12, 12)  # 更大的偏移範圍
        
        # 添加二次趨勢系數
        quadratic_coef = random.uniform(-0.002, 0.002)
        linear_coef = random.uniform(-0.3, 0.3)
        
        # 生成更複雜的信號 
        signal = []
        for t in time_range:
            # 基本正弦波
            sine_component = amplitude * np.sin(frequency * t + phase)
            # 二次趨勢
            trend_component = quadratic_coef * t**2 + linear_coef * t
            # 合併
            value = offset + sine_component + trend_component
            # 添加隨機噪聲
            noisy_value = value + random.uniform(-0.8, 0.8)
            signal.append((t, noisy_value))
        
        # 偶爾添加突變
        if random.random() > 0.7:  # 30% 的機率有突變
            # 在隨機位置添加短暫的突變
            mutation_start = random.randint(num_points // 4, 3 * num_points // 4)
            mutation_length = random.randint(3, 8)
            mutation_amplitude = random.uniform(2, 5) * (1 if random.random() > 0.5 else -1)
            
            for j in range(mutation_start, min(mutation_start + mutation_length, num_points)):
                t, y = signal[j]
                signal[j] = (t, y + mutation_amplitude)
        
        series_list.append(signal)
    
    return series_list

# 生成或載入時間序列數據
def get_time_series(num_series=5, num_points=100, force_new=False):
    data_file = '/Users/lynn/Desktop/time_series_data.pkl'
    
    # 如果檔案已存在且不需要強制重新生成，則載入現有數據
    if os.path.exists(data_file) and not force_new:
        try:
            with open(data_file, 'rb') as f:
                lines = pickle.load(f)
            print(f"載入既有時間序列數據 (共{len(lines)}條曲線)")
            return lines
        except Exception as e:
            print(f"載入數據失敗：{e}，將重新生成")
    
    # 生成新數據
    print(f"生成新的時間序列數據 (共{num_series}條曲線)")
    # 使用不同的隨機種子
    random.seed(123)  # 改變種子值以生成不同的數據
    lines = generate_time_series(num_series, num_points)
    
    # 保存數據以便後續使用
    with open(data_file, 'wb') as f:
        pickle.dump(lines, f)
    
    return lines

# 保存或載入交點信息
def get_intersection_points(lines, force_new=False):
    points_file = '/Users/lynn/Desktop/intersection_points.pkl'
    
    # 如果檔案已存在且不需要強制重新生成，則載入現有數據
    if os.path.exists(points_file) and not force_new:
        try:
            with open(points_file, 'rb') as f:
                data = pickle.load(f)
                all_intersections = data['intersections']
                intersection_info = data['info']
                unified_points = data.get('unified_points', [])
            print(f"載入既有交點數據 (共{len(all_intersections)}個交點)")
            return all_intersections, intersection_info, unified_points
        except Exception as e:
            print(f"載入交點數據失敗：{e}，將重新計算")
    
    # 重新計算所有交點
    print("計算所有交點...")
    all_intersections, intersection_info = find_intersections(lines)
    
    # 保存數據以便後續使用
    with open(points_file, 'wb') as f:
        pickle.dump({'intersections': all_intersections, 'info': intersection_info}, f)
    
    return all_intersections, intersection_info, []

# 允許手動指定特定交點的分類狀態
def manual_classify_points(unified_points, manual_settings=None):
    """
    此函數現已不需要手動分類設定，因為所有先前手動分類的點
    現在都通過identify_good_cross_points函數中的特殊案例分析自動處理。
    保留此函數僅為向後兼容性。
    """
    # 函數現在不做任何事，只返回原始點
    return unified_points

# 主程式開始執行
# 初始化 - 載入或生成數據
lines = get_time_series(force_new=True)  # 設置為True強制重新生成數據

# 找出或載入所有交點
all_intersections, intersection_info, loaded_unified_points = get_intersection_points(lines, force_new=True)
print(f"找到總共 {len(all_intersections)} 個交點")

# 不使用已保存的分類結果，強制重新分類所有交點
unified_points = []
# 識別好的交點，使用更寬鬆的參數設定
good_cross_points, rejected_points = identify_good_cross_points(intersection_info, lines, 
                                                            window_size=3, 
                                                            min_slope_diff=0.2, 
                                                            max_oscillations=4,
                                                            min_correlation_abs=0.3,
                                                            min_slope_significance=0.1)  # 設定最小斜率負值

# 創建一個統一的交點索引，包含所有交點及其分類狀態
point_to_index = {}  # 用於快速查找交點索引

# 首先收集所有良好交點
for point, slope1, slope2, line1_idx, line2_idx, osc1, osc2, corr1, corr2, reason in good_cross_points:
    point_key = f"{point[0]:.4f},{point[1]:.4f}"
    idx = len(unified_points) + 1  # 從1開始的索引
    unified_points.append({
        "index": idx,
        "point": point,
        "status": "good",
        "slopes": (slope1, slope2),
        "slope_diff": abs(slope1 - slope2),
        "oscillations": (osc1, osc2),
        "correlations": (corr1, corr2),
        "lines": (line1_idx, line2_idx),
        "reason": reason
    })
    point_to_index[point_key] = idx

# 然後收集所有被拒絕的交點
for point, slope1, slope2, line1_idx, line2_idx, osc1, osc2, corr1, corr2, reason in rejected_points:
    point_key = f"{point[0]:.4f},{point[1]:.4f}"
    if point_key not in point_to_index:  # 確保沒有重複
        idx = len(unified_points) + 1
        unified_points.append({
            "index": idx,
            "point": point,
            "status": "rejected",
            "slopes": (slope1, slope2),
            "slope_diff": abs(slope1 - slope2),
            "oscillations": (osc1, osc2),
            "correlations": (corr1, corr2),
            "lines": (line1_idx, line2_idx),
            "reason": reason
        })
        point_to_index[point_key] = idx

# 檢查是否還有未分類的交點（這種情況不應該出現，但為安全起見）
for point in all_intersections:
    point_key = f"{point[0]:.4f},{point[1]:.4f}"
    if point_key not in point_to_index:
        idx = len(unified_points) + 1
        unified_points.append({
            "index": idx,
            "point": point,
            "status": "unknown",
            "slopes": (0, 0),
            "slope_diff": 0,
            "oscillations": (0, 0),
            "correlations": (0, 0),
            "lines": (0, 0),
            "reason": "未分類"
        })
        point_to_index[point_key] = idx

# 保存分類結果
with open('/Users/lynn/Desktop/intersection_points.pkl', 'wb') as f:
    pickle.dump({
        'intersections': all_intersections, 
        'info': intersection_info,
        'unified_points': unified_points
    }, f)

# 應用手動分類（此函數現在不再改變任何點的分類）
unified_points = manual_classify_points(unified_points)

print(f"其中 {len([p for p in unified_points if p['status'] == 'good'])} 個為良好交點")
print(f"{len([p for p in unified_points if p['status'] == 'rejected'])} 個點被拒絕")

# 顯示良好交點的詳細資訊
print("\n良好交點詳細資訊:")
print("編號\t時間點\t\t線1斜率\t線2斜率\t斜率差異\t線1震盪\t線2震盪\t線1相關\t線2相關\t線1\t線2")
for p in [p for p in unified_points if p['status'] == 'good']:
    print(f"{p['index']}.\t({p['point'][0]:.2f}, {p['point'][1]:.2f})\t{p['slopes'][0]:.3f}\t{p['slopes'][1]:.3f}\t{p['slope_diff']:.3f}\t{p['oscillations'][0]}\t{p['oscillations'][1]}\t{p['correlations'][0]:.3f}\t{p['correlations'][1]:.3f}\t{p['lines'][0]+1}\t{p['lines'][1]+1}")

# 顯示被拒絕點的摘要
print("\n被拒絕交點的原因摘要:")
rejection_reasons = {}
for p in [p for p in unified_points if p['status'] == 'rejected']:
    reason = p['reason']
    if reason in rejection_reasons:
        rejection_reasons[reason] += 1
    else:
        rejection_reasons[reason] = 1

for reason, count in rejection_reasons.items():
    print(f"{reason}: {count}個點")

# 繪圖
plt.figure(figsize=(12, 8))
plt.grid(True, linestyle='--', alpha=0.7)

# 繪製時間序列
for i, line in enumerate(lines):
    x, y = zip(*line)
    plt.plot(x, y, label=f"Signal {i+1}")

# 繪製所有交點（灰色）
if all_intersections:
    ix, iy = zip(*[p['point'] for p in unified_points])
    plt.scatter(ix, iy, color='gray', s=20, alpha=0.4, label="All Intersections")

# 繪製被拒絕的交點（黃色）
rejected_points_list = [p for p in unified_points if p['status'] == 'rejected']
if rejected_points_list:
    rx, ry = zip(*[p['point'] for p in rejected_points_list])
    plt.scatter(rx, ry, color='orange', s=30, alpha=0.5, label="Rejected Points")

# 繪製好的交點（紅色突出顯示）
good_points_list = [p for p in unified_points if p['status'] == 'good']
if good_points_list:
    gx, gy = zip(*[p['point'] for p in good_points_list])
    plt.scatter(gx, gy, color='red', s=50, label="Good Cross Points")

# 標注所有交點編號
for p in unified_points:
    # 根據狀態選擇顏色
    if p['status'] == "good":
        color = 'red'
        fontweight = 'bold'
    elif p['status'] == "rejected":
        color = 'orange'
        fontweight = 'normal'
    else:
        color = 'gray'
        fontweight = 'normal'
    
    plt.annotate(str(p['index']), (p['point'][0], p['point'][1]), xytext=(5, 5), 
                 textcoords='offset points', color=color, fontweight=fontweight)

plt.title('Time Series Signals with All Intersection Points Labeled')
plt.xlabel('Time (s)')
plt.ylabel('Signal Strength')
plt.legend()

# 設置適當的軸範圍
plt.xlim(0, 100)
# 修復 zip 對象不可索引的問題
y_min = min([min([point[1] for point in line]) for line in lines]) - 1
y_max = max([max([point[1] for point in line]) for line in lines]) + 1
plt.ylim(y_min, y_max)

# 添加標題和網格
plt.tight_layout()
plt.savefig('/Users/lynn/Desktop/all_labeled_intersections.png', dpi=300)
# plt.show()  # 註釋掉，不顯示圖形，避免被中斷

# 輸出所有交點信息
print("\n所有交點信息:")
print("編號\t狀態\t\t時間點\t\t線段信息")

for p in unified_points:
    status_text = "良好" if p['status'] == 'good' else "拒絕" if p['status'] == 'rejected' else "未知"
    print(f"{p['index']}.\t{status_text}\t({p['point'][0]:.2f}, {p['point'][1]:.2f})\t" + 
          f"斜率1: {p['slopes'][0]:.3f}, 斜率2: {p['slopes'][1]:.3f}, " + 
          f"差異: {p['slope_diff']:.3f}, 震盪: {p['oscillations'][0]}/{p['oscillations'][1]}, " + 
          f"相關: {p['correlations'][0]:.3f}/{p['correlations'][1]:.3f}, " + 
          f"線條: {p['lines'][0]+1}/{p['lines'][1]+1}, 原因: {p['reason']}")

# 繪製局部放大圖，展示交點附近的斜率
if good_points_list:
    for p in good_points_list[:min(5, len(good_points_list))]:
        # 創建新圖形
        plt.figure(figsize=(10, 6))
        
        # 計算交點附近的範圍
        x_center = p['point'][0]
        y_center = p['point'][1]
        x_range = 10  # 時間軸的範圍
        
        # 繪製交點附近的線段
        for line_idx, line in enumerate(lines):
            if line_idx == p['lines'][0] or line_idx == p['lines'][1]:  # 只繪製相交的兩條線
                x, y = zip(*line)
                plt.plot(x, y, label=f"Signal {line_idx+1}")
                
                # 標記用於計算斜率的點
                nearby_points = [(j, point) for j, point in enumerate(line) if abs(point[0] - x_center) < x_range]
                if nearby_points:
                    nearby_x = [point[0] for _, point in nearby_points]
                    nearby_y = [point[1] for _, point in nearby_points]
                    plt.scatter(nearby_x, nearby_y, s=30, alpha=0.5)
        
        # 繪製交點
        plt.scatter([x_center], [y_center], color='red', s=100, zorder=5, label=f"Cross Point {p['index']}")
        
        # 添加斜率和震盪信息及相關係數
        plt.text(x_center-x_range*0.4, y_center + 3, f"Slope 1: {p['slopes'][0]:.3f}, Oscil: {p['oscillations'][0]}, Corr: {p['correlations'][0]:.3f}", fontsize=10)
        plt.text(x_center-x_range*0.4, y_center + 2, f"Slope 2: {p['slopes'][1]:.3f}, Oscil: {p['oscillations'][1]}, Corr: {p['correlations'][1]:.3f}", fontsize=10)
        plt.text(x_center-x_range*0.4, y_center + 1, f"Slope Diff: {p['slope_diff']:.3f}", fontsize=10)
        
        # 設置圖形範圍，聚焦於交點附近
        plt.xlim(x_center - x_range, x_center + x_range)
        
        # 動態調整y軸範圍
        nearby_points_all = []
        for line_idx in [p['lines'][0], p['lines'][1]]:
            nearby_points_all.extend([point[1] for j, point in enumerate(lines[line_idx]) if abs(point[0] - x_center) < x_range])
        
        if nearby_points_all:
            y_min_local = min(nearby_points_all) - 1
            y_max_local = max(nearby_points_all) + 3  # 多留空間顯示文字
            plt.ylim(y_min_local, y_max_local)
        
        plt.title(f'Zoomed View of Cross Point {p["index"]}')
        plt.xlabel('Time (s)')
        plt.ylabel('Signal Strength')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'/Users/lynn/Desktop/cross_point_{p["index"]}_detail.png', dpi=300)

print("\n所有圖形已保存至桌面:")
print("- 所有標記交點的圖形: /Users/lynn/Desktop/all_labeled_intersections.png")
for p in good_points_list[:min(5, len(good_points_list))]:
    print(f"- 交點 {p['index']} 的詳細視圖: /Users/lynn/Desktop/cross_point_{p['index']}_detail.png")

# plt.show()  # 註釋掉，不顯示圖形，避免被中斷
