# This is a sample Python script for performing optical flow with intestinal myometer images.
import numpy as np
import cv2 as cv
from alive_progress import alive_bar
import time

file_path = '/Volumes/Backup 2/intestinal-myometer/AMR001--20240417/WIN_20240417_12_51_48_Pro.mp4'
output_path = '/Volumes/Backup 2/intestinal-myometer/AMR001--20240417/WIN_20240417_12_51_48_Pro_tracked.avi'

file_path = '/Volumes/Backup 2/intestinal-myometer/AMR001--20240417/WIN_20240417_12_26_45_Pro_trimmed.mp4'
output_path = '/Volumes/Backup 2/intestinal-myometer/AMR001--20240417/WIN_20240417_12_26_45_Pro_trimmed_tracked.avi'
cap = cv.VideoCapture(file_path)

# params for ShiTomasi corner detection
feature_params = dict(maxCorners=100,
                      qualityLevel=0.3,
                      minDistance=7,
                      blockSize=7)

# Parameters for lucas kanade optical flow
lk_params = dict(winSize=(15, 15),
                 maxLevel=2,
                 criteria=(cv.TERM_CRITERIA_EPS | cv.TERM_CRITERIA_COUNT, 10, 0.03))

# Create some random colors
color = np.random.randint(0, 255, (100, 3))

# Take first frame and find corners in it
ret, old_frame = cap.read()
roi = cv.selectROI("Select the area and press [Enter]", old_frame)

old_frame = old_frame[int(roi[1]):int(roi[1]+roi[3]),
                      int(roi[0]):int(roi[0]+roi[2])]
size = (old_frame.shape[1], old_frame.shape[0])

print("RoI selected")

old_gray = cv.cvtColor(old_frame, cv.COLOR_BGR2GRAY)
p0 = cv.goodFeaturesToTrack(old_gray, mask=None, **feature_params)

# Write optical flow results to file
result = cv.VideoWriter(output_path, cv.VideoWriter_fourcc('I', '4', '2', '0'), 30, size)

# Create a mask image for drawing purposes
mask = np.zeros_like(old_frame)

# Create progress bar
property_id = int(cv.CAP_PROP_FRAME_COUNT)
num_frames = int(cv.VideoCapture.get(cap, property_id))

with alive_bar(num_frames) as bar:
    while True:
        ret, frame = cap.read()
        if not ret:
            break

        frame = frame[int(roi[1]): int(roi[1] + roi[3]), int(roi[0]): int(roi[0] + roi[2])]

        frame_gray = cv.cvtColor(frame, cv.COLOR_BGR2GRAY)

        # calculate optical flow
        p1, st, err = cv.calcOpticalFlowPyrLK(old_gray, frame_gray, p0, None, **lk_params)

        # Select good points
        if p1 is not None:
            good_new = p1[st == 1]
            good_old = p0[st == 1]

        # draw the tracks
        for i, (new, old) in enumerate(zip(good_new, good_old)):
            a, b = new.ravel()
            c, d = old.ravel()
            mask = cv.line(mask, (int(a), int(b)), (int(c), int(d)), color[i].tolist(), 2)
            frame = cv.circle(frame, (int(a), int(b)), 5, color[i].tolist(), -1)
        img = cv.add(frame, mask)

        result.write(img)

        # Now update the previous frame and previous points
        old_gray = frame_gray.copy()
        p0 = good_new.reshape(-1, 1, 2)

        time.sleep(.001)
        bar()

# Release the video read and write objects
result.release()
cap.release()

# Closes all the frames
cv.destroyAllWindows()

print("The video was successfully saved")

print(mask)