import math
import numpy as np
import quaternion
from typing import List

class MadgwickFusion:
    sampling_rate: 0.0
    gyro_error = 0.0 # in rad/s
    beta = 0.0
    zeta = 0.0

    last_acc = List[float]
    last_gyro = List[float]
    last_magneto = List[float]

    estimated_orientation = np.quaternion(0, 0, 0, 0)
    b_x = 0.0 # reference direction of flux in earth frame
    b_z = 0.0 # reference direction of flux in earth frame
    est_gyro_bias = List[float]

    def __init__(self, sampling_rate: float, gyro_err: float, gyro_drift: float):
        self.sampling_rate = sampling_rate
        self.gyro_error = gyro_err

        sqrt_const = math.sqrt(3.0 / 4.0)
        self.beta = gyro_err * sqrt_const
        self.zeta = gyro_drift * sqrt_const

    def filter_update(self, acc, gyro, magneto):
        # normalize data
        acc_norm = self.normalize_vector(acc)
        gyro_norm = self.normalize_vector(gyro)
        magneto = self.normalize_vector(magneto)

    def normalize_vector(self, vector):
        norm = np.linalg.norm(vector, ord=2)
        return [vector[0] / norm, vector[1] / norm, vector[2] / norm]
