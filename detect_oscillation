import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mode
from sklearn.cluster import KMeans
from scipy import stats

def detect_oscillation_regions(signal, tolerance=4.5, min_length=5, debug=False, signal_name=None, trend_threshold=3, strictness_level=1, allowed_consecutive=5, equal_tolerance=0.25):
    """
    Detect oscillation regions in a signal
    
    Parameters:
    - signal: 1D array, signal values
    - tolerance: float, tolerance for direction change
    - min_length: int, minimum length of the interval
    - debug: bool, whether to output debug information
    - signal_name: string, signal name (for debugging)
    - equal_tolerance: float, tolerance for judging value equality
    
    Returns:
    - target_regions: list, containing all detected target intervals [(start, end, interval value list), ...]
    """
    if len(signal) == 0:
        return []
    
    # Convert to numpy array and process only non-zero values
    signal = np.array(signal)
    valid_indices = np.where(signal > 0)[0]
    
    if len(valid_indices) == 0:
        return []
    
    valid_signal = signal[valid_indices]
    
    if debug and signal_name:
        print(f"\n===== {signal_name} Detailed Debug Info =====")
        print(f"Valid data points: {len(valid_indices)}")
        print(f"Unique values: {len(np.unique(valid_signal))}")
        
    # Calculate direction changes between points
    directions = []
    for i in range(1, len(valid_signal)):
        if valid_signal[i] > valid_signal[i-1] + equal_tolerance:
            directions.append("up")
        elif valid_signal[i] < valid_signal[i-1] - equal_tolerance:
            directions.append("down")
        else:
            directions.append("flat")
            
    # Identify upward and downward trend intervals
    trend_regions = []
    current_direction = None
    start_idx = 0
    
    for i, direction in enumerate(directions):
        # To align indices with original valid_indices, add 1
        actual_idx = i + 1
        
        # Start a new trend interval
        if current_direction is None and (direction == "up" or direction == "down"):
            current_direction = direction
            start_idx = i
        # Direction change
        elif current_direction and direction != "flat" and direction != current_direction:
            # Check interval length
            if actual_idx - start_idx >= min_length:
                region_start = valid_indices[start_idx]
                region_end = valid_indices[actual_idx - 1]
                trend_type = "Upward Trend" if current_direction == "up" else "Downward Trend"
                trend_regions.append((region_start, region_end, trend_type))
                
                if debug:
                    print(f"Detected {trend_type}: from time point {region_start}(a) to {region_end}(b), length: {region_end - region_start + 1}")
            
            # Start new trend interval
            current_direction = direction
            start_idx = i
    
    # Process the last trend interval
    if current_direction and len(valid_signal) - start_idx >= min_length:
        region_start = valid_indices[start_idx]
        region_end = valid_indices[-1]
        trend_type = "Upward Trend" if current_direction == "up" else "Downward Trend"
        trend_regions.append((region_start, region_end, trend_type))
        
        if debug:
            print(f"Detected {trend_type}: from time point {region_start}(a) to {region_end}(b), length: {region_end - region_start + 1}")
    
    # Find target intervals (intervals not belonging to trends) based on trend intervals
    target_regions = []
    
    # Sort trend intervals by start time
    trend_regions.sort(key=lambda x: x[0])
    
    # If there are trend intervals, check gaps between trend intervals
    if trend_regions:
        # Check gap before the first trend interval
        if trend_regions[0][0] > valid_indices[0] and trend_regions[0][0] - valid_indices[0] + 1 >= min_length:
            target_start = valid_indices[0]
            target_end = trend_regions[0][0] - 1
            if debug:
                print(f"Detected target region: from time point {target_start}(s) to {target_end}(t), length: {target_end - target_start + 1}")
            target_values = signal[target_start:target_end+1]
            target_values = target_values[target_values > 0]
            target_regions.append((target_start, target_end, target_values))
            
        # Check gaps between trend intervals
        for i in range(len(trend_regions) - 1):
            current_end = trend_regions[i][1]
            next_start = trend_regions[i+1][0]
            
            if next_start > current_end + 1 and next_start - current_end - 1 >= min_length:
                target_start = current_end + 1
                target_end = next_start - 1
                if debug:
                    print(f"Detected target region: from time point {target_start}(s) to {target_end}(t), length: {target_end - target_start + 1}")
                target_values = signal[target_start:target_end+1]
                target_values = target_values[target_values > 0]
                target_regions.append((target_start, target_end, target_values))
        
        # Check gap after the last trend interval
        if trend_regions[-1][1] < valid_indices[-1] and valid_indices[-1] - trend_regions[-1][1] >= min_length:
            target_start = trend_regions[-1][1] + 1
            target_end = valid_indices[-1]
            if debug:
                print(f"Detected target region: from time point {target_start}(s) to {target_end}(t), length: {target_end - target_start + 1}")
            target_values = signal[target_start:target_end+1]
            target_values = target_values[target_values > 0]
            target_regions.append((target_start, target_end, target_values))
    
    # If there are no trend intervals, check if the entire signal can be considered as a target interval
    elif len(valid_indices) >= min_length:
        target_start = valid_indices[0]
        target_end = valid_indices[-1]
        if debug:
            print(f"Detected target region (no trend regions): from time point {target_start}(s) to {target_end}(t), length: {target_end - target_start + 1}")
        target_values = valid_signal
        target_regions.append((target_start, target_end, target_values))
    
    if debug and signal_name:
        print(f"\nNumber of target regions detected: {len(target_regions)}")
        for i, (start, end, values) in enumerate(target_regions):
            print(f"  Region {i+1}: time points {start}(s)-{end}(t) (length: {end-start+1})")
    
    return target_regions

def visualize_oscillation_regions(df, save_path=None, highlight_red_box=False, min_length=5):
    """
    Visualize signal and mark target intervals and trend intervals
    
    Parameters:
    - df: DataFrame, containing signal data
    - save_path: string, path to save image
    - highlight_red_box: bool, whether to highlight red box area
    - min_length: int, minimum length of the interval
    """
    # Modify chart size to make chart taller, not as flat
    plt.figure(figsize=(16, 14), facecolor='whitesmoke')
    ax = plt.gca()
    ax.set_facecolor('whitesmoke')
    
    # Draw grid
    plt.grid(True, linestyle='-', alpha=0.3, color='gray')
    
    # Draw vertical reference lines
    for i in range(0, df.shape[0], 10):
        plt.axvline(x=i, color='darkgray', linestyle='-', alpha=0.4)
    
    # Use bright colors
    colors = ['red', 'blue', 'green', 'purple', 'orange']
    
    # Process each signal, detect and mark target intervals and trend intervals
    for i, col in enumerate(df.columns[1:]):
        signal = df[col].values
        valid_indices = np.where(signal > 0)[0]
        
        if len(valid_indices) == 0:
            continue
        
        color = colors[i % len(colors)]
        
        # Draw full signal
        plt.plot(df['x'][valid_indices], signal[valid_indices], '-o', 
                color=color, markersize=6, 
                label=f"{col}", 
                alpha=0.8, linewidth=1.5)
        
        # Detect target intervals and trend intervals
        debug = True  # Enable debugging to check results
        valid_signal = signal[valid_indices]
        
        # Calculate direction changes between points
        directions = []
        for j in range(1, len(valid_signal)):
            if valid_signal[j] > valid_signal[j-1] + 0.25:
                directions.append("up")
            elif valid_signal[j] < valid_signal[j-1] - 0.25:
                directions.append("down")
            else:
                directions.append("flat")
        
        # Identify trend intervals (upward or downward)
        trend_regions = []
        current_direction = None
        start_idx = 0
        
        for j, direction in enumerate(directions):
            actual_idx = j + 1
            
            # Start a new trend interval
            if current_direction is None and (direction == "up" or direction == "down"):
                current_direction = direction
                start_idx = j
            # Direction change
            elif current_direction and direction != "flat" and direction != current_direction:
                # Check interval length
                if actual_idx - start_idx >= min_length:
                    region_start = valid_indices[start_idx]
                    region_end = valid_indices[actual_idx - 1]
                    trend_type = "Upward Trend" if current_direction == "up" else "Downward Trend"
                    trend_regions.append((region_start, region_end, trend_type))
                
                # Start new trend interval
                current_direction = direction
                start_idx = j
        
        # Process the last trend interval
        if current_direction and len(valid_signal) - start_idx >= min_length:
            region_start = valid_indices[start_idx]
            region_end = valid_indices[-1]
            trend_type = "Upward Trend" if current_direction == "up" else "Downward Trend"
            trend_regions.append((region_start, region_end, trend_type))
            
        # Mark trend interval start and end points a,b (using thin lines)
        for start, end, trend_type in trend_regions:
            # Mark trend start point (a) and end point (b)
            # plt.axvline(x=start, color=color, linestyle=':', alpha=0.8, linewidth=1.5)
            # plt.axvline(x=end, color=color, linestyle=':', alpha=0.8, linewidth=1.5)
            
            # 添加趨勢類型說明（置於a和b之間的上方）
            # mid_x = (start + end) / 2
            # plt.text(mid_x, y_pos_trend, trend_type, color=color, fontsize=10, 
            #         ha='center', bbox=dict(facecolor='white', alpha=0.5, boxstyle='round,pad=0.2'))
            pass
        
        # 獲取目標區間（不屬於趨勢的區間）
        target_regions = detect_oscillation_regions(
            signal, 
            min_length=min_length, 
            debug=debug, 
            signal_name=col
        )
        
        # 標記目標區間的起終點s,t（使用較粗虛線）
        for region_start, region_end, _ in target_regions:
            # 標記起始點(s)和結束點(t)
            plt.axvline(x=region_start, color=color, linestyle='--', alpha=0.8, linewidth=2)
            plt.axvline(x=region_end, color=color, linestyle='--', alpha=0.8, linewidth=2)
            
            # 在虛線上方添加標記「s」(起點)和「t」(終點)
            y_pos = max(signal) + 3 + i * 2  # 錯開不同訊號的標籤高度
            plt.text(region_start, y_pos, f"s", color=color, fontsize=12, fontweight='bold',
                    ha='center', bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.2'))
            plt.text(region_end, y_pos, f"t", color=color, fontsize=12, fontweight='bold',
                    ha='center', bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.2'))
            
            # 添加目標區間說明（置於s和t之間的上方）
            # mid_x = (region_start + region_end) / 2
            # plt.text(mid_x, y_pos, "Target Region", color=color, fontsize=10, 
            #         ha='center', bbox=dict(facecolor='white', alpha=0.5, boxstyle='round,pad=0.2'))
    
    # 設置標題和軸標籤
    plt.title(f'Signal Visualization - Target Region Markers (s,t=Start/End points, Length≥{min_length})', fontsize=16)
    plt.xlabel('Time Points', fontsize=12)
    plt.ylabel('Signal Values', fontsize=12)
    
    # 設定y軸範圍，確保有足夠空間顯示標籤
    y_max = df.iloc[:, 1:].max().max() * 1.3
    plt.ylim(0, max(y_max, 45))
    
    # 添加圖例
    plt.legend(loc='upper right')
    
    # 調整布局
    plt.tight_layout()
    
    # 保存圖形
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Target region boundary image saved to {save_path}")
    
    # 顯示圖形
    plt.show()

# 主函式
if __name__ == "__main__":
    # 讀取數據
    try:
        # 優先使用處理過的數據
        df = pd.read_csv('refined_test_data.csv')
        print("Using processed data: refined_test_data.csv")
    except FileNotFoundError:
        try:
            # 如果找不到處理過的數據，則使用原始測試數據
            df = pd.read_csv('test_data.csv')
            print("Using original test data: test_data.csv")
        except FileNotFoundError:
            print("Data file not found. Please ensure refined_test_data.csv or test_data.csv exists")
            exit(1)
    
    # 定義目標區間檢測參數
    min_length = 5  # 目標區間的最小長度
    
    # 可視化並標記目標區間
    visualize_oscillation_regions(
        df, 
        save_path='signal_visualization.png', 
        highlight_red_box=False,
        min_length=min_length
    ) 
