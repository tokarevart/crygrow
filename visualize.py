import cv2
import numpy as np


def parse_img_data(filename):
    with open(filename, 'r') as file:
        size_line = file.readline()
        l_size = int(size_line.split()[1])

        l_data = []
        for line in file:
            spl = line.split()
            l_data.append(((int(spl[0]), int(spl[1])), (int(spl[2]), int(spl[3]), int(spl[4]))))

        return l_size, l_data


size, data = parse_img_data("automata-image-data.txt")
# BGR palette
img = np.zeros((size, size, 3), np.uint8)
white = (255, 255, 255)
for x in range(size):
    for y in range(size):
        img[x][y] = white

for pack in data:
    img[pack[0][0], pack[0][1]] = (pack[1][2], pack[1][1], pack[1][0])

cv2.imshow("image", img)
cv2.imwrite("fig.png", img)
cv2.waitKey()
