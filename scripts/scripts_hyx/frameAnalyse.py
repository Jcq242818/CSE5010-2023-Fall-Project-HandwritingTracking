import cv2
import numpy as np
import matplotlib.pyplot as plt


def extract_added_black_points(video_path):
    cap = cv2.VideoCapture(video_path)

    # 读取第一帧
    ret, prev_frame = cap.read()
    prev_gray = cv2.cvtColor(prev_frame, cv2.COLOR_BGR2GRAY)

    added_black_points_sequence = []

    while True:
        # 读取下一帧
        ret, frame = cap.read()
        if not ret:
            break

        # 转换为灰度图
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        # 计算相邻帧的差异
        diff = cv2.absdiff(prev_gray, gray)

        # 通过阈值化得到黑点
        _, thresh = cv2.threshold(diff, 30, 255, cv2.THRESH_BINARY)

        # 查找轮廓
        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        # 提取新增黑点的坐标
        added_points = []
        for contour in contours:
            area = cv2.contourArea(contour)
            if area > 5:  # 考虑面积大于5的轮廓为新增黑点
                M = cv2.moments(contour)
                cx = int(M['m10'] / M['m00'])
                cy = int(M['m01'] / M['m00'])
                added_points.append((cx, cy))

        added_black_points_sequence.append(added_points)

        # 更新上一帧
        prev_gray = gray

    cap.release()
    return added_black_points_sequence


def plot_points_sequence(points_sequence):
    fig, ax = plt.subplots()

    for points in points_sequence:
        if points:  # 检查是否为空
            x, y = zip(*points)
            ax.scatter(x, y, marker='.', color='black')  # 使用黑色点表示

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Added Black Points Over Frames')

    plt.show()

def main():
    video_path = '/Users/hyx020222/Documents/GitHub/Passive-Handwriting-Tracking/data/RPReplay_Final1704290401.mov'
    sequence = extract_added_black_points(video_path)
    print(sequence)

    # 删除空时刻
    sequence = [points for points in sequence if points]
    print(sequence)

    plot_points_sequence(sequence)


if __name__ == '__main__':
    main()
