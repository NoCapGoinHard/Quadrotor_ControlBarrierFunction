import numpy as np
from scipy.signal import lti
import numpy as np
from scipy.signal import lsim
import matplotlib.pyplot as plt
from MPCpy.MPC_ObstacleAvoidance import simulation_MPC
from MPCpy.config import *
from MPCpy.MPC_ObstacleFree import simulation_MPC_safe
from math import *
import matplotlib.animation as animation

# set simulation parameters
N_sim = 150
N_mpc = 15
T_s = 0.1
x0 = [8,0,0,0]
delta = 1.5

# select controlled inputs
m = 2 

# define double integrator dynamics (and lti system for simulation)
A, B, C, D = get_sys_with_dim(m)

# define reference vectors
X_ref, Y_ref = get_shifted_circle_ref_trajectory(N_sim, N_mpc, T_s, 4, 0, 4)
U_ref = get_zero_input_ref(N_sim, N_mpc, m)

# define the TOTAL reference trajectory (for all coordinates)
XY_ref = np.vstack((X_ref, Y_ref))

# define time vector (to plot simulation)
t = np.arange(0, (N_sim+N_mpc)*T_s, T_s)

# define obstacle trajectory (moving at crossroads intersection)
X_obs, Y_obs = get_shifted_circle_obs_trajectory(N_sim, N_mpc, T_s, -4, 0, 4)
XY_obs = np.vstack((X_obs, Y_obs))

# simulate tracking with NO obstacle on trajectory
X_sim_NoObs, U_simNoObs = simulation_MPC_safe(A, B, XY_ref, U_ref, N_sim, N_mpc, T_s, x0)

# simulate tracking for current trajectory (with obstacle)
X_sim, U_sim = simulation_MPC(A, B, XY_ref, U_ref, XY_obs, N_sim, N_mpc, T_s, x0, delta)

fig, ax = plt.subplots()

xdata, ydata = [], []
xobs, yobs = [], []
xref, yref = [], []

ln, = ax.plot([], [], '^k')
ln2, = ax.plot([],[], 'b')
ln3, = ax.plot([],[], '--r')

def init():
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    return ln, ln2, ln3

def update(i):
    xdata.append(X_sim[0,i])
    ydata.append(X_sim[2,i])
    '''xdata.append(X_sim_NoObs[0,i])
    ydata.append(X_sim_NoObs[2,i])'''
    xobs.append(XY_obs[0,i])
    yobs.append(XY_obs[1,i])
    xref.append(XY_ref[0,i])
    yref.append(XY_ref[2,i])


    ln.set_data(xdata, ydata)
    ln2.set_data(xobs, yobs)
    ln3.set_data(xref, yref)

    return ln, ln2, ln3

anim = animation.FuncAnimation(fig, update, init_func = init, frames = N_sim, blit = True)

plt.show()