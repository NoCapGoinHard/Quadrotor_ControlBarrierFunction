import numpy as np
import casadi as cs

'''
    Solver CVXPY is not used here because we cannot guarantee that the
    control barrier function is convex (cvxpy only solves convex problems)
'''

'''
    Solver IPOPT is used instead, since it is a more general optimization tool
'''

''' 
    Control barrier function h(X) is defined directly in the MPC problem constraints
    
    ---> h(X[k], X_obs[k], Y[k], Y_obs[k], Z[k], Z_obs[k]) = sqrt( (X[k] - X_obs[k])^2 + (Y[k] - Y_obs[k])^2 + (Z[k] - Z_obs[k])^2 ) - delta >= 0
'''

'''
    The CBF conditon to be applied is (1 dim.):

    ---> h(X[k+1]) - h(X[k]) >= -alpha * h(X[k])

    ---> X[k+1] - X_obs[k+1] - X[k] + X_obs[k] >= -alpha * h(X[k], X_obs[k])
'''

def convex_mpc_quadrotor_safe(A, B, Q, R, X_ref, U_ref, X_pred, x0, N_mpc, dt, iter):

    print("At iteration " + str(iter))
    
    # state and input size
    n, m = B.shape[0], B.shape[1]  

    # define optimizer class and select IPOPT solver
    opt = cs.Opti()
    p_opts, s_opts = {"ipopt.print_level": 0, "expand": True}, {}
    opt.solver("ipopt", p_opts, s_opts)
    
    # decision variables
    X = opt.variable(n, N_mpc)
    U = opt.variable(m, N_mpc-1)

    # objective function (to be minimized)
    objective = 0
    for i in range(N_mpc-1):
        objective += Q*cs.sumsqr(X[:, i]-X_ref[:, i]) + R*cs.sumsqr(U[:, i]-U_ref[:, i])
        
    objective += Q*cs.sumsqr(X[:, N_mpc-1]-X_ref[:, N_mpc-1])

    # constraints (initial condition)
    opt.subject_to( X[:, 0] == x0 )                            
    
    # constraints (dynamics)
    for i in range(N_mpc-1):
        opt.subject_to( X[:, i+1] == X[:,i] + dt * (np.matmul(A, X[:,i]) + np.matmul(B, U[:,i])) )

    # instruct solver to minimze objective function
    opt.minimize(objective)

    # set initial guess to optimize solver
    opt.set_initial(X, X_pred)

    # take solution and return first input value 
    sol = opt.solve()

    return sol.value(U[:,0])


def simulation_MPC_safe(A, B, X_ref, U_ref, sim_time, pred_horizon, sampling_time, initial_condition):
    """Simulation with MPC controller"""

    # Get data
    x0 = initial_condition
    Q =  10
    R = 1
    N_mpc = pred_horizon
    N_sim = sim_time
    dt = sampling_time
    n = 4
    m = 2

    # Simulation
    X_sim = np.zeros((n, N_sim))
    X_sim[:, 0] = x0
    U_sim = np.zeros((m, N_sim-1))

    for i in range(N_sim-1):
        
        # Reference trajectory for current window of N_mpc timesteps
        X_ref_tilde = X_ref[:, i:(i+N_mpc)]
        U_ref_tilde = U_ref[:, i:(i+N_mpc-1)]

        # create custom vector to pass as inital guess to IPOPT
        X_sim_guess = np.zeros((n, N_mpc))
        for j in range(n):
            X_sim_guess[j,0] = X_sim[j, i]

        # Compute optimal input
        U_sim[:, i] = convex_mpc_quadrotor_safe(A, B, Q, R, X_ref_tilde,U_ref_tilde, X_sim_guess, X_sim[:, i],N_mpc, dt, i)

        # Simulate one step
        X_sim[:, i+1] = X_sim[:,i] + dt * (np.matmul(A, X_sim[:,i]) + np.matmul(B, U_sim[:,i]))

        print(str(X_sim[0,i]) + ', ' + str(X_sim[1,i]))
        print(str(X_sim[2,i]) + ', ' + str(X_sim[3,i]))

    return X_sim, U_sim