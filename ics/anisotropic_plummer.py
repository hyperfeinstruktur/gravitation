import numpy as np
from matplotlib import pyplot as plt
import pickle as pkl



####### Model Parameters

G = 4.30091e-3 # pc / solar mass * (kms)^2
a = 10         # pc
M = 6.0e4      # solar mass
q = 1.0
N = 200        # -

####### Generate Positions

def r_of_m(m,a,M):
    return a * ((M/m)**(2./3.) - 1)**(-1./2.)

m_rand = M*np.random.uniform(0.0,1.0,N)
r_rand = r_of_m(m_rand,a,M)
phi_rand = np.random.uniform(0.0,2*np.pi,N)
theta_rand = np.arccos( np.random.uniform(-1.,1.,N) )

x = r_rand * np.sin(theta_rand) * np.cos(phi_rand)
y = r_rand * np.sin(theta_rand) * np.sin(phi_rand)
z = r_rand * np.cos(theta_rand)

X = np.array([x,y,z]).transpose()
print(np.amax(np.sqrt(np.sum(X**2,1))))

####### Generate Velocities
assert q != 0.0, "For q=0, use isotropic_plummer.py"

v_x = 0.0
v_y = 0.0
v_z = 0.0

V = np.array([v_x,v_y,v_z]).transpose()
print(np.amax(V))

####### Generate masses & Pickle results

m = M/N * np.ones(N)
if 1:
    with open('X.pkl','wb') as f:
        pkl.dump(X[:],f)

    with open('V.pkl','wb') as f:
        pkl.dump(V[:],f)

    with open('M.pkl','wb') as f:
        pkl.dump(m[:],f)

####### Optional: Display xy projection
if 0:
    plt.scatter(x,y,s=0.4)
    plt.show()