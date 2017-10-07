import numpy as np
import pandas as pd
import quaternion
from madgwick import MadgwickFusion

def main():
    df = pd.read_csv('sample.csv', header=None, sep=',')
    madgwick = MadgwickFusion(1, 1, 1)

    for row in df.values:
        madgwick.filter_update(row[0:3], row[3:6], row[6:9])


if __name__ == "__main__":
    main()