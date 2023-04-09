import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
import MathOptInterface as MOI
import Ipopt 

import FiniteDiff
import ForwardDiff
import Convex as cvx 
import ECOS
using LinearAlgebra
using Plots
using Random
using JLD2
using Test
import MeshCat as mc 
using Statistics 


include(joinpath(@__DIR__, "utils","fmincon.jl"))


function hermite_simpson(params::NamedTuple, x1::Vector, x2::Vector, u, dt::Real)::Vector
    """ We use the discretized version of the dynamics to create a dynamics constraint  """
    
    x_half = ( (1 / 2) * (x1 + x2) ) + ( (dt / 8 ) * ( combined_dynamics(params, x1, u) - combined_dynamics(params, x2, u) ) )
    
    ẋ1 = combined_dynamics(params, x1, u)
    ẋ2 = combined_dynamics(params, x2, u)
    ẋ_half = combined_dynamics(params, x_half, u)
    
    resid = x1 + ( (dt/6) * (ẋ1 + ( 4 * ẋ_half ) + ẋ2)) - x2
    return resid

end



function quadrotor_dynamics_constraints(params::NamedTuple, Z::Vector)::Vector
    """ Use the dynamics to create the dynamics constraints """

    idx, N, dt = params.idx, params.N, params.dt
        
    c = zeros(eltype(Z), idx.nc)
    
    for i = 1:(N-1)
        
        #this is the state for all three quadcopters
        xi = Z[idx.x[i]]
        
        #this is the control for all three quadcopters
        ui = Z[idx.u[i]] 
        
        xi1p = Z[idx.x[i+1]]
        
        # discretized dynamics
        c[idx.c[i]] = hermite_simpson(params, xi, xi1p, ui, dt)
    end
    return c 
end



function quadrotor_cost(params::NamedTuple, Z::Vector)::Real
    idx, N = params.idx, params.N
    x1g, x2g, x3g = params.x1g, params.x2g, params.x3g
    Q, R, Qf = params.Q, params.R, params.Qf
    
    # stack the goals 
    xg = [x1g; x2g; x3g]
    
    J = 0 
    for i = 1:(N-1)
        xi = Z[idx.x[i]]
        ui = Z[idx.u[i]]
        
        # stage cost
        J += 0.5 * transpose(xi - xg) * Q * (xi - xg) + 0.5 * transpose(ui) * R * ui

    end
    
    # dont forget terminal cost 
    xn = Z[idx.x[N]]
    
    J += 0.5 * transpose(xn - xg) * Qf * (xn - xg) 
    
    return J 
end


function quadrotor_equality_constraint(params::NamedTuple, Z::Vector)::Vector
    N, idx = params.N, params.idx
    x1ic, x2ic, x3ic = params.x1ic, params.x2ic, params.x3ic
    x1g, x2g, x3g = params.x1g, params.x2g, params.x3g

    # TODO: return all of the equality constraints 
    
    # initial condition constraints
    ic1_ceq = Z[idx.x[1]][1:6] - x1ic 
    ic2_ceq = Z[idx.x[1]][7:12] - x2ic
    ic3_ceq = Z[idx.x[1]][13:18] - x3ic
    
    ic_ceq = [ic1_ceq; ic2_ceq; ic3_ceq]
    
    # goal condition constraints
    goal1_ceq = Z[idx.x[N]][1:6] - params.x1g
    goal2_ceq = Z[idx.x[N]][7:12] - params.x2g
    goal3_ceq = Z[idx.x[N]][13:18] - params.x3g

    goal_ceq = [goal1_ceq; goal2_ceq; goal3_ceq]
    
    # dynamics constraint
    dyn_ceq = quadrotor_dynamics_constraints(params, Z)
    
    # stack them all
    c_eq = [ic_ceq ; goal_ceq; dyn_ceq]
    
    return c_eq 
end


function quadrotor_ineq_constraint(params, Z)
    
    idx, N = params.idx, params.N 
    
    # collision constraints
    quad1pos = [Z[idx.x[i]][1:2] for i = 1:N]
    quad2pos = [Z[idx.x[i]][7:8] for i = 1:N]
    quad3pos = [Z[idx.x[i]][13:14] for i = 1:N]

    quad1_quad2 = norm.(quad1pos .- quad2pos).^2
    quad2_quad3 = norm.(quad2pos .- quad3pos).^2
    quad1_quad3 = norm.(quad1pos .- quad3pos).^2

    c_ineq = [quad1_quad2; quad1_quad3; quad2_quad3]

    return c_ineq
end


function create_idx(nx,nu,N)
    # This function creates some useful indexing tools for Z 
    # x_i = Z[idx.x[i]]
    # u_i = Z[idx.u[i]]
    
    # Feel free to use/not use anything here.
    # our Z vector is [x0, u0, x1, u1, …, xN]
    nz = (N-1) * nu + N * nx # length of Z 
    x = [(i - 1) * (nx + nu) .+ (1 : nx) for i = 1:N]
    u = [(i - 1) * (nx + nu) .+ ((nx + 1):(nx + nu)) for i = 1:(N - 1)]
    
    # constraint indexing for the (N-1) dynamics constraints when stacked up
    c = [(i - 1) * (nx) .+ (1 : nx) for i = 1:(N - 1)]
    nc = (N - 1) * nx # (N-1)*nx 
    
    return (nx=nx,nu=nu,N=N,nz=nz,nc=nc,x= x,u = u,c = c)
end

"""
    quadrotor_reorient

Function for returning collision free trajectories for 3 quadrotors. 

Outputs:
    x1::Vector{Vector}  # state trajectory for quad 1 
    x2::Vector{Vector}  # state trajectory for quad 2 
    x3::Vector{Vector}  # state trajectory for quad 3 
    u1::Vector{Vector}  # control trajectory for quad 1 
    u2::Vector{Vector}  # control trajectory for quad 2 
    u3::Vector{Vector}  # control trajectory for quad 3 
    t_vec::Vector
    params::NamedTuple

The resulting trajectories should have dt=0.2, tf = 5.0, N = 26
where all the x's are length 26, and the u's are length 25. 

Each trajectory for quad k should start at `xkic`, and should finish near 
`xkg`. The distances between each quad should be greater than 0.8 meters at 
every knot point in the trajectory. 
"""
function quadrotor_reorient(;verbose=true)
    
    # problem size 
    nx = 18 
    nu = 6
    dt = 0.2
    tf = 5.0 
    t_vec = 0:dt:tf 
    N = length(t_vec)
    
    # indexing 
    idx = create_idx(nx,nu,N)
    
    # initial conditions and goal states 
    lo = 0.5 
    mid = 2 
    hi = 3.5 
    x1ic = [-2,lo,0,0,0,0]  # ic for quad 1 
    x2ic = [-2,mid,0,0,0,0] # ic for quad 2 
    x3ic = [-2,hi,0,0,0,0]  # ic for quad 3 

    x1g = [2,mid,0,0,0,0]   # goal for quad 1 
    x2g = [2,hi,0,0,0,0]    # goal for quad 2 
    x3g = [2,lo,0,0,0,0]    # goal for quad 3 
    
    # LQR cost 
    Q = diagm(ones(nx))
    R = 0.1*diagm(ones(nu))
    Qf = 10*diagm(ones(nx))
    
    # load all useful things into params 
    # TODO: include anything you would need for a cost function (like a Q, R, Qf if you were doing an 
    # LQR cost)
    params = (x1ic=x1ic,
              x2ic=x2ic,
              x3ic=x3ic,
              x1g = x1g,
              x2g = x2g,
              x3g = x3g,
              dt = dt,
              N = N,
              idx = idx,
              Q = Q,
              R = R,
              Qf = Qf,
              mass = 1.0, # quadrotor mass 
              g = 9.81,   # gravity 
              ℓ = 0.3,    # quadrotor length 
              J = .018)   # quadrotor moment of inertia 
    
    # TODO: solve for the three collision free trajectories however you like
    
    # create a reference for all paths
    ref = [range(x1ic, x1g, length=N), range(x2ic, x2g, length=N), range(x3ic, x3g, length=N)]
    x0 = 0.001 * randn(idx.nz)
    
    # replace the states with the warm-start ref 
    for i = 1:N
        x0[idx.x[i]] = [ref[1][i]; ref[2][i]; ref[3][i]]
    end
    
    # no bounds on the primal variable
    x_l = -Inf * ones(idx.nz)
    x_u = Inf * ones(idx.nz)
    
    # since I created the bounds directly the lower bound is zero
    c_l = 0.8^2 * ones(3*N)
    c_u = Inf * ones(3*N)
    
    # choose diff type (try :auto, then use :finite if :auto doesn't work)
    diff_type = :auto 

    Z = fmincon(quadrotor_cost,quadrotor_equality_constraint, quadrotor_ineq_constraint,x_l,x_u,c_l,c_u,x0,
    params, diff_type; tol = 1e-6, c_tol = 1e-6, max_iters = 10_000, verbose = verbose)
    
    # return the trajectories 
    
    x1 = [ Z[idx.x[i]][1:6] for i = 1:N]
    x2 = [ Z[idx.x[i]][7:12]  for i = 1:N]
    x3 = [ Z[idx.x[i]][13:18] for i = 1:N]
                
    u1 = [ Z[idx.u[i]][1:2] for i = 1:(N-1)]
    u2 = [ Z[idx.u[i]][3:4] for i = 1:(N-1)]
    u3 = [ Z[idx.u[i]][5:6] for i = 1:(N-1)]
                            
    return x1, x2, x3, u1, u2, u3, t_vec, params 
end
    