'''This module collects some extra functions that are used by the test
drivers (the other files of the form "test_*.py").'''

import autograd.numpy as np

def get_linear_speed_tau(vx, vy):
    s = lambda x, y: 1.0/(1.0 + vx*x + vy*y)
    return lambda x, y: np.arccosh(
        1 + s(x, y)*(vx**2 + vy**2)*(x**2 + y**2)/2
    )/np.sqrt(vx**2 + vy**2)

def get_linear_speed_s(vx, vy):
    return lambda x, y: 1.0/(1.0 + vx*x + vy*y)
