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

'''
    GOOD BEHAVIOUR at crossroads: N_mpc=5, delta = 0.5 (Q=R=1)

    

    !!! DATA for T_impact = 2 (alpha = pi/8) and T_s = 0.1 !!!
    ACCEPTABLE BEHAVIOUR at circle: M_mpc=5, delta=0.25 (Q=R=1) (NO CBF)

    GOOD BEHAVIOUR at circle: N_mpc=15, delta=0.5 (but limit N_sim to 38) (Q=R=1) (NO CBF)

    GOOD BEHAVIOUR at circle: N_mpc=10, delta=1 (Q=R=1) (NO CBF)

    GOOD BEHAVIOUR at circle: N_mpc=10, delta=1 / 0.9 (Q=R=1) (CBF)

    VERY GOOD BEHAVIOUR at circle: N_mpc=15, delta=1.5 (Q=R=1) (CBF) (even better for Q=10, R=1)

    !!! DATA for T_impact = 5 (alpha = pi/10) and T_s = 0.1 !!!
    INTERESTING BEHAVIOUR at N_mpc = 5, delta = 0.75 (Q=10, R=1) (CBF)

    VERY GOOD BEHAVIOUR at N_mpc=5, delta = 0.75 (Q=R=1) (CBF)

    VERY COOL BEHAVIOUR AT N_mpc=15, delta=0.2 (Q=R=1) (CBF)

    VERY GOOD BEHAVIOUR at N_mpc=15, delta = 1 (Q=R=1) (CBF)

    INTERESTING BEHAVIOUR at N_mpc=15, delta=1 (Q=10, R=1) (CBF)

'''

# set simulation parameters
N_sim = 150
N_mpc = 10
T_s = 0.1
x0 = [0,0,0,0]
delta = 0.2

# select controlled inputs
m = 2 

# define double integrator dynamics (and lti system for simulation)
A, B, C, D = get_sys_with_dim(m)

# define reference vectors
X_ref = get_x_ref_trajectory(N_sim, N_mpc, T_s)
Y_ref = get_y_ref_trajectory(N_sim, N_mpc, T_s)
U_ref = get_zero_input_ref(N_sim, N_mpc, m)

# define the TOTAL reference trajectory (for all coordinates)
XY_ref = np.vstack((X_ref, Y_ref))

# define time vector (to plot simulation)
t = np.arange(0, (N_sim+N_mpc)*T_s, T_s)

# define obstacle trajectory (moving at crossroads intersection)
X_obs, Y_obs = get_circle_obs_trajectory(N_sim, N_mpc, T_s, 5)
XY_obs = np.vstack((X_obs, Y_obs))

# simulate tracking with NO obstacle on trajectory
X_sim_NoObs, U_simNoObs = simulation_MPC_safe(A, B, XY_ref, U_ref, N_sim, N_mpc, T_s, x0)

# simulate tracking for current trajectory (with obstacle)
X_sim, U_sim = simulation_MPC(A, B, XY_ref, U_ref, XY_obs, N_sim, N_mpc, T_s, x0, delta)

fig, ax = plt.subplots()

xdata, ydata = [], []
xobs, yobs = [], []
xref, yref = [], []

ln, = ax.plot([], [], 'r', label='MPC')
ln2, = ax.plot([],[], 'b', label='Obstacle')
ln3, = ax.plot([],[], '--r', label='Reference')

leg = ax.legend(loc="lower left")

def init():
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 8)
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