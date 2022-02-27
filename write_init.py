import numpy as np
import pickle as pkl

# ==== Notes
# Units:
# distance: kiloparsec
# mass: solar mass
# velocity: km/s
# (G = 4.30091e-6 # pc / Ms * (kms)^2)
# Time: seconds

inputfname = 'big_plummer.in'
# Parameters as in the init file
year_in_sec = 86400
tFin=9e9 * year_in_sec # 1 Gyr
G= 6.674000e-11 # 4.30091e-6        # 
nb_obj=4000
nsteps=2000
output='big_plummer.out'
sampling=5
adaptatif=0
epsilon=1.000000e-03


# Nx3 arrays
x = np.empty((nb_obj,3))
v = np.empty((nb_obj,3))
m = np.empty(nb_obj)

# x = np.array([[1,1,2],[0,0,1]])
# v = np.array([[10,6,1],[8,6,2]])
with open('X.pkl','rb') as f:
    x = pkl.load(f)
with open('V.pkl','rb') as f:
    v = pkl.load(f)
with open('M.pkl','rb') as f:
    m = pkl.load(f)

f = open(inputfname,'w')



# Box bounds
inf= -400. #-7.714250e+17
sup= 400. #7.714250e+17

def write_to_init(parname,val):
    f.write(parname+'='+str(val) + '\n')

write_to_init('tFin',tFin)
write_to_init('G',G)
write_to_init('nb_obj',nb_obj)
write_to_init('nsteps',nsteps)
write_to_init('output',output)
write_to_init('sampling',sampling)
write_to_init('adaptatif',adaptatif)
write_to_init('epsilon',epsilon)
write_to_init('inf',inf)
write_to_init('sup',sup)

# Write positions
for i in range(len(m)):
    y = x[i,:]
    s = 'y' + str(i+1)
    write_to_init(s + '_1',y[0])
    write_to_init(s + '_2',y[1])
    write_to_init(s + '_3',y[2])

# Write velocities
for i in range(len(m)):
    vel = v[i,:]
    s = 'v' + str(i+1)
    write_to_init(s + '_1',vel[0])
    write_to_init(s + '_2',vel[1])
    write_to_init(s + '_3',vel[2])

# Write masses
for i in range(len(m)):
    s = 'm' + str(i+1)
    write_to_init(s,m[i])


f.close()