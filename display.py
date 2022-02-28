import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
#import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

outpath = 'big_plummer.out'
out = np.loadtxt(outpath)
nb_obj = 2000
nsteps = 2000
sampling = 5
nframes = int(nsteps / sampling)
lim = 30
#print(out)
print(nframes)

t = out[:,-1] / 86400. / 1e9
init = out[0,:]
#print(len(init))
xi = init[0:3*nb_obj:3]
yi = init[1:3*nb_obj:3]
zi = init[2:3*nb_obj:3]

print(xi.shape)
#print(len(xi))
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(projection='3d',fc='black',box_aspect=(1,1,1))
sc = ax.scatter3D(xi,yi,zi,s=1,c='white')
ax.set_xlim(-lim,lim)
ax.set_ylim(-lim,lim)
ax.set_zlim(-lim,lim)
ax.set_title('$t=0$',fontsize=25)
ax.spines['bottom'].set_color('white')
ax.spines['top'].set_color('white') 
ax.spines['right'].set_color('white')
ax.spines['left'].set_color('white')
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')
ax.tick_params(axis='z', colors='white')
ax.set_xlabel('$x$ [pc]',color='white')
ax.set_ylabel('$y$ [pc]',color='white')
ax.set_zlabel('$z$ [pc]',color='white')

def update(k):
    #print(k)
    inst = out[k,:]
    #print(inst[0])
    #print(np.shape(inst))
    pos = np.array([inst[0:3*nb_obj:3],inst[1:3*nb_obj:3],inst[2:3*nb_obj:3]])
    #sc.set_offsets(pos)
    sc._offsets3d = (inst[0:3*nb_obj:3],inst[1:3*nb_obj:3],inst[2:3*nb_obj:3])
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.set_zlim(-lim,lim)
    ax.set_title('$t=$'+format(t[k],'.1f') + ' Gyrs',fontsize=25)


ani = FuncAnimation(fig, update,frames=nframes, interval=30, repeat=1) 

plt.show()