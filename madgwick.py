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

    SE_q = np.quaternion(0, 0, 0, 0)
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
        magneto_norm = self.normalize_vector(magneto)

        # compute the objective function and jacobian
        SE_q_comp = self.SE_q.components
        f = np.array(
            [2 * SE_q_comp[0] * SE_q_comp[3] - 2 * SE_q_comp[0] * SE_q_comp[2] - acc_norm[0],
             2 * SE_q_comp[1] * SE_q_comp[1] + 2 * SE_q_comp[2] * SE_q_comp[3] - acc_norm[1],
             1 - 2 * SE_q_comp[1] * SE_q_comp[1] - 2 * SE_q_comp[2] * SE_q_comp[2] - acc_norm[2],
             2 * self.b_x * (0.5 - SE_q_comp[2] * SE_q_comp[2] - SE_q_comp[3] * SE_q_comp[3]) + 2 * self.b_z * (SE_q_comp[1] * SE_q_comp[3] - SE_q_comp[0] * SE_q_comp[2]) - magneto_norm[0],
             2 * self.b_x * (SE_q_comp[1] * SE_q_comp[2] - SE_q_comp[0] * SE_q_comp[3] + 2 * self.b_z * (SE_q_comp[0] * SE_q_comp[1] + SE_q_comp[2] * SE_q_comp[3]) - magneto_norm[1]),
             2 * self.b_x * (SE_q_comp[0] * SE_q_comp[3]) + 2 * self.b_z * (0.5 - SE_q_comp[1] * SE_q_comp[1] - SE_q_comp[2] * SE_q_comp[2]) - magneto_norm[2]])

    def normalize_vector(self, vector):
        norm = np.linalg.norm(vector, ord=2)
        return [vector[0] / norm, vector[1] / norm, vector[2] / norm]
