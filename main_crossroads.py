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
'''

# set simulation parameters
N_sim = 250
N_mpc = 5
T_s = 0.1
x0 = [0,0,0,0]
delta = 1

# select controlled inputs
m = 2 

# define double integrator dynamics (and lti system for simulation)
A, B, C, D = get_sys_with_dim(m)

# define reference vectors
X_ref = get_x_ref_trajectory(N_sim, N_mpc, T_s)
Y_ref = get_y_ref_trajectory(N_sim, N_mpc, T_s)
U_ref = get_zero_input_ref(N_sim, N_mpc, m)

# define time vector (to plot simulation)
t = np.arange(0, (N_sim+N_mpc)*T_s, T_s)

# define obstacle trajectory (fixed)
'''X_obs = get_fixed_obstacle_traj(N_sim, N_mpc, 0)
Y_obs = get_fixed_obstacle_traj(N_sim, N_mpc, 5)
XY_obs = np.vstack((X_obs, Y_obs))'''

# define obstacle trajectory (moving at crossroads intersection)
X_obs, Y_obs = get_crossroads_obs_traj(N_sim, N_mpc, T_s, 5)
XY_obs = np.vstack((X_obs, Y_obs))

# define the TOTAL reference trajectory (for all coordinates)
XY_ref = np.vstack((X_ref, Y_ref))

# simulate tracking with NO obstacle on trajectory
X_sim_NoObs, U_simNoObs = simulation_MPC_safe(A, B, XY_ref, U_ref, N_sim, N_mpc, T_s, x0)

# create a reference (nominal) input for obstacle case
'''U_ref_non_zero = get_u_ref_trajectory(N_sim, N_mpc, m, U_simNoObs)'''

# simulate tracking for current trajectory
X_sim, U_sim = simulation_MPC(A, B, XY_ref, U_ref, XY_obs, N_sim, N_mpc, T_s, x0, delta)
'''X_sim, U_sim = simulation_MPC(A, B, XY_ref, U_ref_non_zero, XY_obs, N_sim, N_mpc, T_s, x0, delta)'''

'''# plot results (with time on x axis)
plt.plot(t[0:N_sim], X_sim[0,:], 'b')
plt.plot(t[0:N_sim], X_ref[0,0:N_sim], 'r--')

plt.plot(t[0:N_sim], X_sim[1,:], 'k')
plt.plot(t[0:N_sim], X_ref[1,0:N_sim], 'g--')'''

# plot results for 2 dimensions
'''plt.plot(X_sim[0,:], X_sim[2,:], 'b')

plt.grid(alpha=0.3)
plt.xlabel('x')
plt.show()'''



fig, ax = plt.subplots()
xdata, ydata = [], []
xobs, yobs = [], []
ln, = ax.plot([], [], 'r')
ln2, = ax.plot([],[], 'b')

def init():
    ax.set_xlim(-0.1, 0.1)
    ax.set_ylim(0, 16)
    return ln,

def update(i):
    xdata.append(X_sim[0,i])
    ydata.append(X_sim[2,i])
    xobs.append(XY_obs[0,i])
    yobs.append(XY_obs[1,i])

    ln.set_data(xdata, ydata)
    ln2.set_data(xobs, yobs)

    return ln, ln2

anim = animation.FuncAnimation(fig, update, init_func = init, frames = N_sim, blit = True)

plt.show()
    

