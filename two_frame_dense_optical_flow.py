import cv2 as cv
import numpy as np
from os import listdir
from os.path import isfile, join
import subprocess
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
import os
import glob

def plot_quiver(ax, flow, spacing, margin=0, **kwargs):
    """Plots less dense quiver field.

    Args:
        ax: Matplotlib axis
        flow: motion vectors
        spacing: space (px) between each arrow in grid
        margin: width (px) of enclosing region without arrows
        kwargs: quiver kwargs (default: angles="xy", scale_units="xy")
    """
    h, w, *_ = flow.shape

    nx = int((w - 2 * margin) / spacing)
    ny = int((h - 2 * margin) / spacing)

    x = np.linspace(margin, w - margin - 1, nx, dtype=np.int64)
    y = np.linspace(margin, h - margin - 1, ny, dtype=np.int64)

    flow = flow[np.ix_(y, x)]
    u = flow[:, :, 0]
    v = flow[:, :, 1]

    kwargs = {**dict(angles="xy", scale_units="xy"), **kwargs}
    ax.quiver(x, y, u, v, **kwargs)

    ax.set_ylim(sorted(ax.get_ylim(), reverse=True))
    ax.set_aspect("equal")


def generate_quiver_video(only_files, output_path, vid_name):
    image_old = cv.imread(join(folder_src, only_files[0]))
    image_old = cv.cvtColor(image_old, cv.COLOR_BGR2GRAY)

    fig, ax = plt.subplots()
    loopnum = 0
    for i in range(1, len(only_files), 50):
        image_new = cv.imread(join(folder_src, only_files[i]))
        image_new = cv.cvtColor(image_new, cv.COLOR_BGR2GRAY)

        flow = cv.calcOpticalFlowFarneback(image_old, image_new, None, 0.5, 3, 15, 3, 5, 1.2, 0)

        plt.imshow(image_new, cmap='gray')
        plot_quiver(ax, flow, spacing=40, scale=1, color="#ff44ff")
        plt.savefig(output_path + "/file%02d.png" % loopnum)
        loopnum += 1

        image_old = image_new

    os.chdir(output_path)
    subprocess.call([
        'ffmpeg', '-framerate', '8', '-i', 'file%02d.png', '-r', '30', '-pix_fmt', 'yuv420p',
        vid_name
    ])

    for file_name in glob.glob("*.png"):
        os.remove(file_name)


folder_main = '/Users/jdow403/Desktop/AWB015_VID006'
folder_src = folder_main + '/images/'
folder_outputs = folder_main+'/outputs'

only_files = [f for f in listdir(folder_src) if (isfile(join(folder_src, f)) and '.Bmp' in f)]

generate_quiver_video(only_files, folder_outputs, 'AWB015_quiver.mp4')
