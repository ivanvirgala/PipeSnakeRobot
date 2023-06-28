import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

tStart = 0
tStop = 60
dt = 0.01

t = np.arange(tStart, tStop, dt)

def model_diff(x, t):
    