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


lambda0 = 1.55e-6              # free space wavelength (m)
c0 = 3e8                    # speed of light in vacuum (m/s)
omega = 2*np.pi*c0/lambda0  # angular frequency (2pi/s)
dl = 0.05                  # grid size (L0)
NPML = [20, 20]             # number of pml grid points on x and y borders
pol = 'Ez'                  # polarization (either 'Hz' or 'Ez')
source_amp = 1e-9           # amplitude of modal source (make around 1 for nonlinear effects)

# define material constants
n_index = 2.8              # refractive index
eps_m = n_index**2          # relative permittivity

# geometric parameters for a 1 -> 2 port device
L = 20     # length of box (L0) = micron 
N = 10         # Num output ports 
H = 20    # height of box (L0)
w = .5        # width of waveguides (L0)
d = H/10     # distance between waveguides (L0)
l = 5         # length of waveguide from PML to box (L0)
spc = 3       # space between box and PML (L0)

# define permittivity of three port system
eps_r, design_region = N_IO_port(N, L, H, w, d, l, spc, dl, NPML, eps_m)
(Nx, Ny) = eps_r.shape
nx, ny = int(Nx/2), int(Ny/2)            # halfway grid points

# make a new simulation object
simulation = Simulation(omega, eps_r, dl, NPML, pol)

print("Computed a domain with {} grids in x and {} grids in y".format(Nx,Ny))
print("The simulation has {} grids per free space wavelength".format(int(lambda0/dl/simulation.L0)))

# plot the permittivity distribution
simulation.plt_eps(outline=False)
plt.show()

# i = 9 # Waveguide num (start at 0)
# y_i =  # Waveguide y-pos
for i in range(N):
    # set the input waveguide modal source
    simulation.add_mode(neff=np.sqrt(eps_m), direction_normal='x', center=[NPML[0]+int(l/2/dl), ny-((float(i)-float(N-1)/2.0)*d/dl + 0.5)], width=int(H/13/dl), scale=source_amp, order=1)
simulation.setup_modes()
print('Setup sources!')

for i in range(N):
    simulation.modes[i].scale = 1e-7
simulation.setup_modes()
print('Set new instance of source amps!')

(Hx, Hy, Ez) = simulation.solve_fields()
print('Solved forward fields!')

for i in range(N):
    simulation.modes[i].scale = 1e-6
simulation.setup_modes()
print('Set new instance of source amps!')

(Hx, Hy, Ez) = simulation.solve_fields()
print('Solved forward fields for instance 2!')

# for i in range(N):
#     simulation.add_mode(neff=np.sqrt(eps_m), direction_normal='x', center=[-NPML[0]-int(l/2/dl), ny-((float(i)-float(N-1)/2.0)*d/dl + 0.5)], width=int(H/13/dl), scale=source_amp)
# simulation.setup_modes()
# print('Solved adjoint fields!')