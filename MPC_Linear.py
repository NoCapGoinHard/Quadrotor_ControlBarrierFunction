import numpy as np
import cvxpy as cp


def convex_mpc_quadrotor(A, B, Q, R, X_ref, U_ref, x0, N_mpc, dt, iter):
    
    n, m = B.shape[0], B.shape[1]  # State and input size

    # Decision variables
    X = cp.Variable((n, N_mpc))
    U = cp.Variable((m, N_mpc-1))

    # Objective function (quadratic)
    objective = 0
    for i in range(N_mpc-1):
        objective += 0.5*cp.quad_form(X[:, i]-X_ref[:, i], Q) + 0.5*cp.quad_form(U[:, i]-U_ref[:, i], R)

    objective += 0.5*cp.quad_form(X[:, N_mpc-1]-X_ref[:, N_mpc-1], Q)

    # Constraints
    constraints = [X[:, 0] == x0]                            # Initial condition
    #constraints += [X[:,N_mpc-1] == X_ref[:,0]]
    for i in range(N_mpc-1):
        constraints += [X[:, i+1] == X[:,i] + dt * (np.matmul(A, X[:,i]) + np.matmul(B, U[:,i]))]  # Dynamics

    prob = cp.Problem(cp.Minimize(objective), constraints)
    prob.solve(verbose=False)
    return U.value[:, 0]


def simulation_MPC(A, B, X_ref, U_ref, sim_time, pred_horizon, sampling_time, initial_condition):
    """Simulation with MPC controller."""

    # Get data
    x0 = initial_condition
    Q =  10*np.identity(2)
    R = 0.1*np.identity(1)
    N_mpc = pred_horizon
    N_sim = sim_time
    dt = sampling_time
    n = 2
    m = 1

    # Simulation
    X_sim = np.zeros((n, N_sim))
    X_sim[:, 0] = x0
    U_sim = np.zeros((m, N_sim-1))

    for i in range(N_sim-1):
        
        # Reference trajectory for current window of N_mpc timesteps
        X_ref_tilde = X_ref[:, i:(i+N_mpc)]
        U_ref_tilde = U_ref[:, i:(i+N_mpc-1)]

        # Compute optimal input
        U_sim[:, i] = convex_mpc_quadrotor(A, B, Q, R, X_ref_tilde, U_ref_tilde, X_sim[:, i],N_mpc, dt, i)

        # Simulate one step
        X_sim[:, i+1] = X_sim[:,i] + dt * (np.matmul(A, X_sim[:,i]) + np.matmul(B, U_sim[:,i]))

    return X_sim, U_sim

