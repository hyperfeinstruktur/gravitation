import numpy as np
from matplotlib import pyplot as plt
import pickle as pkl



####### Model Parameters

G = 4.30091e-3 # pc / solar mass * (kms)^2
a = 10         # pc
M = 6.0e4      # solar mass
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

def relative_potential(r,a,M): # = - phi
    return G*M / np.sqrt(r**2 + a**2)

def relative_energy(r,v,a,M):
    return relative_potential(r,a,M) - 0.5 * v**2

def DF(r,v,a,M):
    A = 24. * np.sqrt(2) / (7*np.pi**3)
    eps = relative_energy(r,v,a,M)
    if eps < 0: return 0
    else : return A * G**(-5.) * M**(-4.) * a**2 * eps**(3.5)

DF_vect = np.vectorize(DF)

def v_max(r,a,M):
    return np.sqrt(2.*relative_potential(r,a,M))

def sample_velocity(r,a,M):
    v_e = v_max(r,a,M)
    # x = 0.
    # y = 0.
    # while True:
    #     x = v_e*np.random.uniform(0.,1.)
    #     y = DF(r,0.,a,M)*np.random.uniform(0.,1.) # The DF is maximal at r=0 (analytic)
    #     if y < DF(r,x,a,M) and x < v_max(r,a,M): return x,y
    while True:
        X4 = np.random.uniform()
        X5 = np.random.uniform()
        if 0.1*X5 < (X4**2*(1.-X4**2)**(3.5)): return X4*v_e, X5*v_e

vel_rand = np.empty(N)
for i in range(len(vel_rand)):
    vel_rand[i] = sample_velocity(r_rand[i],a,M)[0]

phi_rand = np.random.uniform(0.0,2*np.pi,N)
theta_rand = np.arccos( np.random.uniform(-1.,1.,N) )

v_x = vel_rand * np.sin(theta_rand) * np.cos(phi_rand)
v_y = vel_rand * np.sin(theta_rand) * np.sin(phi_rand)
v_z = vel_rand * np.cos(theta_rand)

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
if 1:
    plt.scatter(x,y,s=0.4)
    plt.show()