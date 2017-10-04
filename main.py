import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import quaternion

def main():
    Fuse('sample.csv')

    df = pd.read_csv(file, header=None, sep=',')

    quaternions_x = []
    quaternions_y = []
    quaternions_z = []
    quaternions_w = []

    for i in range(0, df[0].size):
        accel = np.array([float(df[0][i]), float(df[1][i]), float(df[2][i])])
        gyro = np.array([float(df[3][i]), float(df[4][i]), float(df[5][i])])
        magneto = np.array([float(df[6][i]), float(df[7][i]), float(df[8][i])])


        quaternions_x.append(x)
        quaternions_y.append(y)
        quaternions_z.append(z)
        quaternions_w.append(w)


if __name__ == "__main__":
    main()