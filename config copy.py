import numpy as np

def get_one_dim_system():
    A = np.matrix(' 0 1; 0 0')
    B = np.matrix(' 0; 1')
    C = np.matrix('1 0; 0 1')
    D = np.matrix('0; 0')

    return A, B, C, D

def get_two_dim_system():
    A = np.matrix(' 0 1 0 0 ; 0 0 0 0; 0 0 0 1; 0 0 0 0')
    B = np.matrix(' 0 0; 1 0; 0 0; 0 1')
    C = np.matrix('1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1')
    D = np.matrix('0 0; 0 0; 0 0; 0 0')

    return A, B, C, D

def get_sys_with_dim(m):
    if m == 1:
        return get_one_dim_system()
    elif m == 2:
        return get_two_dim_system()
    else:
        exit('Must select dimension bewtween 1 and 2')

def get_fixed_obstacle_traj(N_sim, N_mpc, coord):
    
    '''
    INPUT: 
        N_sim: number of simulations instants
        N_mpc: prediction horizon
        x_coord: coordinate of fixed obstacle
        
    OUTPUT:
        X_obs: vector of size (N_mpc+N_sim) containing trajectory of the obstacle in time'''
    
    X_obs = np.zeros((1, N_sim+N_mpc))

    for i in range(N_sim+N_mpc):
        X_obs[0,i] = coord

    return X_obs

def get_crossroads_obs_traj(N_sim, N_mpc, T_s, y_coord):
    
    # define x coordinate for obstacle trajectory
    X_obs = np.zeros((1, N_sim+N_mpc))

    t = np.arange(0,  (N_sim+N_mpc)*T_s, T_s)  

    traj = t

    for i in range(N_sim+N_mpc):
        X_obs[:,i] = traj[i]

    X_obs = X_obs - y_coord

    # define y coordinate 
    Y_obs = np.zeros((1, N_sim+N_mpc)) + y_coord

    # define the complete obstacle trajectory
    XY_obs = np.hstack((X_obs, Y_obs))

    return X_obs, Y_obs

def get_x_ref_trajectory(N_sim, N_mpc, T_s):
    
    '''
    INPUT: 
        N_sim: number of simulations instants
        N_mpc: prediction horizon
        T_s: sampling time
        
    OUTPUT:
        X_ref: vector of size (N_mpc+N_sim) containing reference trajectory'''


    R = 1;      # radious of the circumpherence
    w = 1;      # angular velocity
    r=0.15      # shady circumf. parameter
    
    t = np.arange(0,  (N_sim+N_mpc)*T_s, T_s)       
    X_ref = np.zeros((2, N_sim+N_mpc))

    '''traj = np.zeros(1, len(t))
    traj_dot  = np.zeros(1, len(t))

    # fill reference vectors with trajectory to track
    for i in range(N_sim+N_mpc):
        X_ref[:,i] = [traj[i], traj_dot[i]]'''

    return X_ref

def get_y_ref_trajectory(N_sim, N_mpc, T_s):
    
    '''
    INPUT: 
        N_sim: number of simulations instants
        N_mpc: prediction horizon
        T_s: sampling time
        
    OUTPUT:
        X_ref: vector of size (N_mpc+N_sim) containing reference trajectory'''
    
    t = np.arange(0,  (N_sim+N_mpc)*T_s, T_s)       
    Y_ref = np.zeros((2, N_sim+N_mpc))

    traj = t
    traj_dot  = np.ones_like(t)

    # fill reference vectors with trajectory to track
    for i in range(N_sim+N_mpc):
        Y_ref[:,i] = [traj[i], traj_dot[i]]

    return Y_ref

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

def get_u_ref_trajectory(N_sim, N_mpc, m, U):
    '''
    INPUT: 
        N_sim: number of simulations instants
        N_mpc: prediction horizon
        m: dimension of the system (2 dim. for our simulations)
        U: vector containig reference (nominal) input to be elongated

    OUTPUT:
        U_ref: vector of size (N_mpc+N_sim) containing reference input'''


    U_ref = np.zeros((m, N_sim+N_mpc))

    for i in range(N_mpc-1):
        for j in range(m):
            U_ref[j,i] = U[j,i]

    return U_ref

def get_zero_input_ref(N_sim, N_mpc, m):

    '''
    INPUT: 
        N_sim: number of simulations instants
        N_mpc: prediction horizon
        
    OUTPUT:
        U_ref: vector of size (N_mpc+N_sim) containing reference input'''

    U_ref = np.zeros((m, N_sim+N_mpc))
   
    return U_ref

