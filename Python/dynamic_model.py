import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

tStart = 0
tStop = 60
dt = 0.01

t = np.arange(tStart, tStop, dt)

# =============================== Variables declaration ===============================
N = 15                          # Number of snake robot links
l = 0.0525                      # 1/2 of snake robot link length
m = 0.406                       # weight of snake robot link
g = 9.81                        # gravitational acceleration
diameter = 0.3                  # diameter of pipeline
d = diameter - 2*l              # diameter for snake robot link
pipelineLength = 3              # length of pipeline
umax = 3                        # maximum value of joint torque
qmax = 400*dt                   # maximum angular velocity of joint
dt = 0.01                       # increment of time
Erub = 400000                   # constant for contact modeling
vrub = 0.49                     # constant for contact modeling
friction = 0                    # Friction: 0 - Coulomb friction, 1 - viscous friction
ct = 0.015                      # viscous friction coefficient in tangential direction
cn = 0.03                       # viscous friction coefficient in normal direction
ctPipe = 0.08                   # viscous friction coefficient of pipeline side walls
ut = 0.15                       # Coulomb friction coefficient in tangetial direction
un = 0.3                        # Coulomb friction coefficient in normal direction
utPipe = 0.2                    # Coulomb friction coefficient of pipeline side walls
contact = 0                     # Side walls: 0 - without  contact, 1 - with contact
minLinkVel = 0.001              # minimum link velocity for contact modeling
kp = 25                         # constant for P controller
kd = 10                         # constant for D controller
alpha = 0.3981                  # reference trajectory parameter
omega = 0.6936                  # reference trajectory parameter
delta = 0.4914                  # reference trajectory parameter

# Auxiliary vectors & matrices
A = np.zeros((N,N+1))
D = np.zeros(())
K = np.zeros((N,N+1))
J = np.zeros((N,N+1))
J2 = np.zeros((N,N-1))
J4 = np.zeros((N,N-1))
e = np.ones((N+1,1))
e[N+1][1] = 0


for i1 in range(N):
    for i2 in range(N+1):
        if(i1==i2):
            J[i1][i2] = -1
        if(i2==i1+1):
            J[i1][i2] = 1

for i1 in range(N):
    for i2 in range(N+1):
        if(i1==i2):
            K[i1][i2] = 1

for i1 in range(N):
    for i2 in range(N-1):
        if(i1==i2):
            J2[i1][i2] = 1
        if(i2==i1+1):
            J2[i1][i2] = 1

for i1 in range(N):
    for i2 in range(N-1):
        if(i1==i2):
            J4[i1][i2] = 1
        if(i2==i1+1):
            J4[i1][i2] = -1

J3 = -J2
J1 = -J4

# =====================================================================================

#def model_diff(x, t):
    