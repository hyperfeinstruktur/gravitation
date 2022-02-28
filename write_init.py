import numpy as np
import pickle as pkl

def write_to_init(file,parname,val):
    file.write(parname+'='+str(val) + '\n')

# == Units ==:
# distance: parsec
# mass: solar mass
# velocity: km/s
# (G = 4.30091e-3 # pc / Ms * (kms)^2)
# Time: seconds

# Name of input file to be generated
inputfname = 'plummer.in'

# Simulation parameters
year_in_sec = 86400
tFin=8e9 * year_in_sec
nb_obj=200
self_gravity = 1

# Numerical Parameters
nsteps=1500
output='plummer.out'
sampling=5


# Box bounds (parsecs)
inf= -500.
sup= 500.

# Initialize Initial condition vectors
x = np.empty((nb_obj,3)) # Nx3
v = np.empty((nb_obj,3)) # Nx3
m = np.empty(nb_obj)     # Nx1

# Import ICs
with open('X.pkl','rb') as f:
    x = pkl.load(f)
with open('V.pkl','rb') as f:
    v = pkl.load(f)
with open('M.pkl','rb') as f:
    m = pkl.load(f)

# Write to input file
f = open(inputfname,'w')
write_to_init(f,'tFin',tFin)
write_to_init(f,'nb_obj',nb_obj)
write_to_init(f,'nsteps',nsteps)
write_to_init(f,'output',output)
write_to_init(f,'sampling',sampling)
write_to_init(f,'inf',inf)
write_to_init(f,'sup',sup)
write_to_init(f,'self_gravity',self_gravity)

# Write positions
for i in range(len(m)):
    y = x[i,:]
    s = 'y' + str(i+1)
    write_to_init(f,s + '_1',y[0])
    write_to_init(f,s + '_2',y[1])
    write_to_init(f,s + '_3',y[2])

# Write velocities
for i in range(len(m)):
    vel = v[i,:]
    s = 'v' + str(i+1)
    write_to_init(f,s + '_1',vel[0])
    write_to_init(f,s + '_2',vel[1])
    write_to_init(f,s + '_3',vel[2])

# Write masses
for i in range(len(m)):
    s = 'm' + str(i+1)
    write_to_init(f,s,m[i])

f.close()