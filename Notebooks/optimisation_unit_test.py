# Main code for optimisation - unit test 
# Author: Joel Sved 
# Date: 3/7/23

# 1. Define import libraries etc
import numpy as np
import matplotlib.pylab as plt
import copy

# add angler to path (not necessary if pip installed)
import sys
sys.path.append("..")

# import the main simulation and optimization classes
from angler import Simulation, Optimization

# import some structure generators
from angler.structures import three_port, two_port, N_port, N_IO_port
from angler.filter import rho2eps, get_W, rho2rhot, rhot2rhob, rhob2eps
import seaborn as sns 

# LaTeX plot formatting
#Matplotlib params for figures
params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,}
plt.rcParams.update(params)
plt.rcParams['text.usetex'] = True
plt.rc('text.latex', preamble=r'\usepackage[cm]{sfmath}')
sns.set_context("talk",font_scale=1)

###########################################################
# 2. Define constants for simulation 
lambda0 = 1.55e-6           # free space wavelength (m)
c0 = 3e8                    # speed of light in vacuum (m/s)
omega = 2*np.pi*c0/lambda0  # angular frequency (2pi/s)
dl = 0.05                   # grid size (L0)
NPML = [20, 20]             # number of pml grid points on x and y borders
pol = 'Ez'                  # polarization (either 'Hz' or 'Ez')
source_amp = 1e-3           # mW? amplitude of modal source (make around 1 for nonlinear effects)
n_index = 2.8               # refractive index (2-D Si)
eps_m = n_index**2          # relative permittivity

# Define geomtric properties for a N -> N port device
L = 20              # length of box (L0) = micron 
N = 10              # Num output ports 
H = 20              # height of box (L0)
w = .5              # width of waveguides (L0)
d = H/10            # distance between waveguides (L0)
l = 5               # length of waveguide from PML to box (L0)
spc = 3             # space between box and PML (L0)
width= int(H/13/dl) # Width of source 
R = 5               # radius convolution

###########################################################
# 3. Configure simulation objects 
eps_r, design_region = N_IO_port(N, L, H, w, d, l, spc, dl, NPML, eps_m)
(Nx, Ny) = eps_r.shape
nx, ny = int(Nx/2), int(Ny/2)            # halfway grid points

# Make new simulation objects, 1 fwd and 1 adj
simulation_fwd = Simulation(omega, eps_r, dl, NPML, pol)
simulation_adj = Simulation(omega, eps_r, dl, NPML, pol)

print("Computed a domain with {} grids in x and {} grids in y".format(Nx,Ny))
print("The simulation has {} grids per free space wavelength".format(int(lambda0/dl/simulation_fwd.L0)))

# Plot the permittivity distribution
simulation_fwd.plt_eps(outline=False)
plt.show()

W = get_W(Nx, Ny, design_region, NPML, R=R) # Called by rho2eps function 

###########################################################
# 4. Configure source objects 
# Setup mode source params
src_pos_fwd = []
src_pos_adj = []

for i in range(N):
    center_fwd = [NPML[0]+int(l/2/dl), ny-((float(i)-float(N-1)/2.0)*d/dl + w)]
    center_adj = [-NPML[0]-int(l/2/dl), ny-((float(i)-float(N-1)/2.0)*d/dl + w)]
    # set the input waveguide modal source
    simulation_fwd.add_mode(neff=np.sqrt(eps_m), direction_normal='x', center=center_fwd, width=width, scale=source_amp)
    simulation_adj.add_mode(neff=np.sqrt(eps_m), direction_normal='x', center=center_adj, width=width, scale=source_amp)
    src_pos_fwd.append(center_fwd)
    src_pos_adj.append(center_adj)

srcs_fwd = simulation_fwd.setup_modes()
srcs_adj = simulation_adj.setup_modes()

###########################################################
# 5. Load dataset class (for now just a single instance of each)
#TODO :: load class 

