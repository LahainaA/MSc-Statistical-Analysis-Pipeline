import numpy as np
import matplotlib.pyplot as plt
import os

# 1. CONFIGURATION
DURATION_SEC = 3.7
TARGET_FPS = 60
TOTAL_FRAMES = int(DURATION_SEC * TARGET_FPS)
PATH_WIDTH_PX = 800
SINE_AMPLITUDE_RATIO = 0.05
SINE_CYCLES = 3
ARC_STEPS = 500

# 2. REPLICATE THE ALGORITHM 
def get_point(p_norm, width, amp_ratio, cycles):
    x = p_norm * width
    amp_px = width * amp_ratio
    angle = cycles * p_norm * 2 * np.pi
    y = -amp_px * np.sin(angle)
    return np.array([x, y])

def build_lut(width, amp_ratio, cycles, steps):
    segments = []
    total_length = 0
    prev_point = get_point(0, width, amp_ratio, cycles)
    for i in range(1, steps + 1):
        p = i / steps
        curr_point = get_point(p, width, amp_ratio, cycles)
        dist = np.linalg.norm(curr_point - prev_point)
        total_length += dist
        segments.append((p, total_length))
        prev_point = curr_point
    return segments, total_length

def get_p_at_distance(target_dist, lut, total_length):
    if target_dist <= 0: return 0.0
    if target_dist >= total_length: return 1.0
    for i, (p_end, dist_end) in enumerate(lut):
        if dist_end >= target_dist:
            if i == 0: p_start, dist_start = 0.0, 0.0
            else: p_start, dist_start = lut[i-1]
            fraction = (target_dist - dist_start) / (dist_end - dist_start)
            return p_start + fraction * (p_end - p_start)
    return 1.0

# 3. SIMULATE 
print("Running validation simulation...")
lut, total_path_length = build_lut(PATH_WIDTH_PX, SINE_AMPLITUDE_RATIO, SINE_CYCLES, ARC_STEPS)
speed_per_frame = total_path_length / TOTAL_FRAMES
frame_coordinates = []
for frame_num in range(TOTAL_FRAMES + 1):
    target_distance = frame_num * speed_per_frame
    p = get_p_at_distance(target_distance, lut, total_path_length)
    frame_coordinates.append(get_point(p, PATH_WIDTH_PX, SINE_AMPLITUDE_RATIO, SINE_CYCLES))

# 4. ANALYZE 
steps = np.diff(np.array(frame_coordinates), axis=0)
step_sizes = np.linalg.norm(steps, axis=1)
mean_step = np.mean(step_sizes)
max_dev = ((np.max(step_sizes) - mean_step) / mean_step) * 100
cv = np.std(step_sizes) / mean_step

print("\n" + "="*40)
print(f"RESULTS FOR THESIS:")
print(f"Max Deviation: {max_dev:.4f}%")
print(f"Coef. Variation: {cv:.5f}")
print("="*40)

#  5. PLOT 
plt.figure(figsize=(10, 5))
plt.plot(step_sizes, label='Step Size')
plt.axhline(mean_step, color='r', linestyle='--', label='Mean')
plt.title(f'Kinematic Validation (Max Dev: {max_dev:.2f}%)')
plt.xlabel('Frame Number')
plt.ylabel('Pixel Step Size')

# Save to Desktop
save_path = os.path.join(os.path.expanduser("~"), "Desktop", "arc_length_validation_plot.png")
plt.savefig(save_path)
print(f"Plot saved to: {save_path}")
