import math
import numpy as np
from scipy import optimize
from pyquaternion import Quaternion

from typing import List

class MadgwickFusion:
    sampling_rate: 0.0
    gyro_error = 0.0 # in rad/s
    beta = 0.0
    zeta = 0.0

    last_acc = List[float]
    last_gyro = List[float]
    last_magneto = List[float]

    SE_q = np.array([0, 0, 0, 0], np.float)
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
        acc_norm = np.insert(self.normalize_vector(acc), 0, 0)
        gyro_norm = np.insert(self.normalize_vector(gyro), 0, 0)
        magneto_norm = np.insert(self.normalize_vector(magneto), 0, 0)

        # Group 1 magnetic distortion compensation recalculation
        E_h = self.SE_q * np.transpose(magneto_norm) * np.transpose(Quaternion(self.SE_q).conjugate.elements)
        E_b = np.array([0, np.linalg.norm(np.array([E_h[0], E_h[1]]), ord=2), 0, E_h[3]])

        # magnetometer and accelerometer fusion

        # # compute SE_q using (14)
        # optimize.minimize(lambda SE_q, E_d, S_s: np.tensordot((np.tensordot(SE_q, E_d, axes=0)), S_s, axes=0), self.SE_q, )
        #
        # f = np.array(
        #     [2 * self.SE_q[0] * self.SE_q[3] - 2 * self.SE_q[0] * self.SE_q[2] - acc_norm[0],
        #      2 * self.SE_q[1] * self.SE_q[1] + 2 * self.SE_q[2] * self.SE_q[3] - acc_norm[1],
        #      1 - 2 * self.SE_q[1] * self.SE_q[1] - 2 * self.SE_q[2] * self.SE_q[2] - acc_norm[2],
        #      2 * self.b_x * (0.5 - self.SE_q[2] * self.SE_q[2] - self.SE_q[3] * self.SE_q[3]) + 2 * self.b_z * (self.SE_q[1] * self.SE_q[3] - self.SE_q[0] * self.SE_q[2]) - magneto_norm[0],
        #      2 * self.b_x * (self.SE_q[1] * self.SE_q[2] - self.SE_q[0] * self.SE_q[3] + 2 * self.b_z * (self.SE_q[0] * self.SE_q[1] + self.SE_q[2] * self.SE_q[3]) - magneto_norm[1]),
        #      2 * self.b_x * (self.SE_q[0] * self.SE_q[3]) + 2 * self.b_z * (0.5 - self.SE_q[1] * self.SE_q[1] - self.SE_q[2] * self.SE_q[2]) - magneto_norm[2]])

    def normalize_vector(self, vector):
        norm = np.linalg.norm(vector, ord=2)
        return [vector[0] / norm, vector[1] / norm, vector[2] / norm]

    def J_g(self, SE_q):
        return np.array([
            [-2 * SE_q[2], 2 * SE_q[3], -2 * SE_q[0], 2 * SE_q[1]],
            [2 * SE_q[1], 2 * SE_q[0], 2 * SE_q[3], 2 * SE_q[2]],
            [0 -4 * SE_q[1], -4 * SE_q[2], 0]
        ])

    def J_b(self, SE_q, E_b):
        return np.array([
            [-2 * E_b[2] * SE_q[2], 2 * E_b[2] * SE_q[3], -4 * E_b[0] * SE_q[2] - 2 * E_b[2] * SE_q[0],
             -4 * E_b[0] * SE_q[3] + 2 * E_b[2] * SE_q[1]],
            [-2 * E_b[0] * SE_q[3] + 2 * E_b[2] * SE_q[1], 2 * E_b[0] * SE_q[2] + 2 * E_b[2] * SE_q[0],
             2 * E_b[0] * SE_q[1] + 2 * E_b[2] * SE_q[3], -2 * E_b[0] * SE_q[0] + 2 * E_b[2] * SE_q[2]],
            [2 * E_b[0] * SE_q[2], 2 * E_b[0] * SE_q[3] - 4 * E_b[2] * SE_q[1], 2 * E_b[0] * SE_q[0] -
             4 * E_b[2] * SE_q[2], 2 * E_b[0] * SE_q[1]]
        ])

    def f_g(self, SE_q, S_a):
        return np.array([
            [2 * (SE_q[1] * SE_q[3] - SE_q[0] * SE_q[2]) - S_a[0]],
            [2 * (SE_q[0] * SE_q[1] + SE_q[2] * SE_q[3]) - S_a[1]],
            [2 * (0.5 - SE_q[1]**2 - SE_q[2]**2) - S_a[2]]
        ])

    def f_b(self, SE_q, E_b, S_m):
        return np.array(
            [[2 * E_b[0] * (0.5 - SE_q[2]**2 - SE_q[3]**2) + 2 * E_b[2] * (SE_q[1] * SE_q[3] - SE_q[0] * SE_q[2]) -
              S_m[0]],
             [2 * E_b[0] * (SE_q[1] * SE_q[2] - SE_q[0]*SE_q[3]) + 2 * E_b[2] * (SE_q[0] * SE_q[1] + SE_q[2] * SE_q[3]) -
              S_m[1]],
             [2 * E_b[0] * (SE_q[0] * SE_q[2] + SE_q[1] * SE_q[3]) + 2 * E_b[2] * (0.5 - SE_q[1]**2 - SE_q[2]**2)
              - S_m[1]]])

    def Jg