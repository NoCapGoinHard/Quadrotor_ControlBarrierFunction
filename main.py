import pandas as pd
import numpy as np
from MPC_KF_Obstacle import simulation_MPC_KF
from MPC_ObstacleAvoidance import simulation_MPC
from config import *
from math import *
import matplotlib.animation as animation
import matplotlib.pyplot as plt

def get_circle_obs_trajectory(N_sim, N_mpc, T_s, t_impact):

    '''
    INPUT: 
        N_sim: number of simulations instants
        N_mpc: prediction horizon
        T_s: sampling time
        t_imapct: instant of time when we want the impact to happen
        
    OUTPUT:
        X_ref: vector of size (N_mpc+N_sim) containing reference cricle trajectory'''
    
    t = np.arange(0, (N_sim+N_mpc)*T_s, T_s)   
    
    # pick circumerence radius
    r = t_impact

    # particulare case for r = 2
    X_obs = r*np.cos((np.pi/(2*t_impact))*t)  
    Y_obs = r*np.sin((np.pi/(2*t_impact))*t)

    return X_obs, Y_obs
def get_shifted_circle_ref_trajectory(N_sim, N_mpc, T_s, x0, y0, r):
    '''
    INPUT: 
        N_sim: number of simulations instants
        N_mpc: prediction horizon
        T_s: sampling time
        x0: amount of x-axis shift
        y0: amount of y-axis shift 
        r: circle radius
        
    OUTPUT:
        X_ref: vector of size (N_mpc+N_sim) containing reference cricle trajectory'''
    
    t = np.arange(0, (N_sim+N_mpc)*T_s, T_s)   

    X_ref = np.zeros((2, N_sim+N_mpc))
    Y_ref = np.zeros((2, N_sim+N_mpc))  

    # define trajectory to track
    x_traj = x0 + r*np.cos(t) 
    x_traj_dot = -r*np.sin(t)

    y_traj = y0 + r*np.sin(t)
    y_traj_dot = r*np.cos(t)

    # fill reference vectors with trajectory to track
    for i in range(N_sim+N_mpc):
        X_ref[:,i] = [x_traj[i], x_traj_dot[i]]

    for i in range(N_sim+N_mpc):
        Y_ref[:,i] = [y_traj[i], y_traj_dot[i]]

    return X_ref, Y_ref
def get_shifted_circle_obs_trajectory(N_sim, N_mpc, T_s, x0, y0, r):
    '''
    INPUT: 
        N_sim: number of simulations instants
        N_mpc: prediction horizon
        T_s: sampling time
        x0: amount of x-axis shift
        y0: amount of y-axis shift 
        r: circle radius
        
    OUTPUT:
        X_ref: vector of size (N_mpc+N_sim) containing reference cricle trajectory'''
    
    t = np.arange(0, (N_sim+N_mpc)*T_s, T_s)   
    
    X_obs = x0 - r*np.cos(t)  
    Y_obs = y0 + r*np.sin(t)

    return X_obs, Y_obs

# read csv file containing predictions
df = pd.read_csv('../MPC_predictions_opp.csv')
df = df.reset_index()

temp = df.to_numpy()


# select useful rows and correspinding values from dataframe
predictions = temp[:,1:len(temp[0])]

# set simulation parameters
N_sim = 150
N_mpc = 5
T_s = 0.1
x0 = [8,0,0,0]
delta = 0.5
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

# get x and y predictions at each timestep
X_obs_pred = np.zeros((N_sim-1, N_mpc+1))
for i in range(N_sim-1):
    X_obs_pred[i,0] = X_obs[i]
    for j in range(N_mpc):
        X_obs_pred[i, j+1] = predictions[i,(j)*3]
X_obs_pred = np.vstack((X_obs_pred[0,:], X_obs_pred))

Y_obs_pred = np.zeros((N_sim-1, N_mpc+1))
for i in range(N_sim-1):
    Y_obs_pred[i,0] = Y_obs[i]
    for j in range(N_mpc):
        Y_obs_pred[i, j+1] = predictions[i,(j)*3+1]
Y_obs_pred = np.vstack((Y_obs_pred[0,:], Y_obs_pred))


# stack efficiently for MPC for loop
XY_obs_pred = np.zeros((2*N_sim, N_mpc+1))
for i in range(N_sim):
    XY_obs_pred[2*i:2*i+2, :] = np.vstack((X_obs_pred[i,:], Y_obs_pred[i,:]))


# simulate tracking with NO obstacle on trajectory
X_sim_nominal, U_sim_nominal = simulation_MPC(A, B, XY_ref, U_ref, XY_obs, N_sim, N_mpc, T_s, x0, delta)

# simulate tracking for current trajectory (with obstacle)
X_sim, U_sim = simulation_MPC_KF(A, B, XY_ref, U_ref, XY_obs_pred, N_sim, N_mpc, T_s, x0, delta)

fig, ax = plt.subplots()

xdata, ydata = [], []
xobs, yobs = [], []
xref, yref = [], []
xpred, ypred = [], []

ln, = ax.plot([], [], 'r', label='MPC+KF')
ln2, = ax.plot([],[], 'b', label='Obstacle')
ln3, = ax.plot([],[], '--r', label='MPC')
ln4, = ax.plot([],[], 'og', label='Predictions')

leg = ax.legend(loc="lower left")

def init():
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    return ln, ln2, ln3

def update(i):
    xdata.append(X_sim[0,i])
    ydata.append(X_sim[2,i])
    xobs.append(XY_obs[0,i])
    yobs.append(XY_obs[1,i])
    xpred=X_obs_pred[i, 0:N_mpc]
    ypred=Y_obs_pred[i, 0:N_mpc]
    xref.append(X_sim_nominal[0,i])
    yref.append(X_sim_nominal[2,i])
    ln.set_data(xdata, ydata)
    ln2.set_data(xobs, yobs)
    ln3.set_data(xref, yref)
    ln4.set_data(xpred, ypred)

    return ln, ln2, ln3, ln4

anim = animation.FuncAnimation(fig, update, init_func = init, frames = N_sim, blit = True)

plt.show()













