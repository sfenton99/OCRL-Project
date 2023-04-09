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



function cost(params::NamedTuple, Z::Vector)::Real
    idx, N = params.idx, params.N
    xg = params.xg
    Q, R, Qf = params.Q, params.R, params.Qf

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


function equality_constraint(params::NamedTuple, Z::Vector)::Vector
    N, idx = params.N, params.idx
    xic = params.xic
    xg = params.xg

    # initial condition constraints
    ic_ceq = Z[idx.x[1]] - xic
    
    # goal condition constraints
    goal_ceq = Z[idx.x[N]] - params.xg

    # dynamics constraint
    dyn_ceq = quadrotor_dynamics_constraints(params, Z)
    
    # stack them all
    c_eq = [ic_ceq ; goal_ceq; dyn_ceq]
    
    return c_eq 
end


function quadrotor_ineq_constraint(params, Z)
    idx, N = params.idx, params.N

    # collision constraints
    pos = [Z[idx.x[i]][1:2] for i = 1:N]

    c_obs1 = norm.(pos .- pos_obs1).^2
    c_obs2 = norm.(pos .- pos_obs2).^2
    c_obs3 = norm.(pos .- pos_obs3).^2

    c_ineq = [c_obs1; c_obs2; c_obs3]

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
function generate_trajectory(;verbose=true)
    
    # problem size 
    nx = 3
    nu = 2
    dt = 0.2
    tf = 10.0
    t_vec = 0:dt:tf 
    N = length(t_vec)
    
    # indexing
    idx = create_idx(nx,nu,N)
    
    # initial conditions and goal states 
    lo = 0.5 
    mid = 2 
    hi = 3.5 
    xic = [-2,lo,0,0,0,0]  # ic
    xg = [2,mid,0,0,0,0]   # goal
    
    # LQR cost 
    Q = diagm(ones(nx))
    R = 0.1*diagm(ones(nu))
    Qf = 10*diagm(ones(nx))
    
    # load all useful things into params 
    params = (xic=xic,
              xg = xg,
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

    # create a reference for all paths
    ref = range(xic, xg, length=N)
    x0 = 0.001 * randn(idx.nz)
    # replace the states with the warm-start ref 
    for i = 1:N
        x0[idx.x[i]] = ref[i]
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
    x = [ Z[idx.x[i]] for i = 1:N]
    u = [ Z[idx.u[i]] for i = 1:(N-1)]
                            
    return x, u, t_vec, params
end
    