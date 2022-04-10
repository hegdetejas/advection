import numpy as np
import matplotlib.pyplot as plt
import ast


def extract_data(filename):
    f = open(filename, "r")

    file_str = f.read()
    data = ast.literal_eval(file_str)
    NT = data[0][0]
    grid = data[1]
    np_grid = np.array(grid)

    return (NT, np_grid)


def visualize(filename, NT, np_grid):
    fig = plt.figure(figsize=(10, 8))

    ax = fig.add_subplot(111)
    if NT == 0:
        title = "NT = 0 (Gaussian Initialisation)"
    else:
        title = "NT = " + str(NT)
    ax.set_title(title, fontsize=20)
    ax.set_ylabel('y (cm)')
    ax.set_xlabel('x (cm)')
    plt.imshow(np_grid)
    ax.set_aspect('equal')
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

    cax = fig.add_axes([0.13, 0.11, 0.85, 0.77])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical', ax=ax)
    plt.savefig(filename)
    plt.clf()


if __name__ == "__main__":
    data = extract_data("../out/out_0.txt")
    visualize("../out/out_0.png", data[0], data[1])

    data = extract_data("../out/out_NTby2.txt")
    visualize("../out/out_NTby2.png", data[0], data[1])

    data = extract_data("../out/out_NT.txt")
    visualize("../out/out_NT.png", data[0], data[1])
