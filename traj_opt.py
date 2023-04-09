
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

model = pyo.ConcreteModel()

def dynamics(params, x, u):
    return xdot

def hermite_simpson(params, x1, x2, u, dt):
    return xnext

def cost(params,Z):
    return J

def dynamics_constraints(params, Z):
    return c

def equality_constraints(params, Z):
    return c

def inequality_constraints(params, Z):
    return c

def create_idx(nx,nu,N):
    return (nx=nx,nu=nu,N=N,nz=nz,nc=nc,x= x,u = u,c = c)

def traj_opt():
    nx =
    nu =
    dt = 0.2
    tf = 10.0
    t_vec =
    N = length(t_vec)
    idx = create_idx(nx, nu, N)



    return x, u, t_vec, params


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
