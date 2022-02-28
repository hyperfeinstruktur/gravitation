import numpy as np
import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.animation import FuncAnimation

# Path to output file to be displayed
outpath = 'plummer.out'

# Simulation Parameters
nb_obj = 200
nsteps = 1500
sampling = 5

# Figure parameters
lim = 30
save_vid = 0
figsize = 8

# Loading initial frame
out = np.loadtxt(outpath)
t = out[:,-1] / 86400. / 1e9
init = out[0,:]
xi = init[0:3*nb_obj:3]
yi = init[1:3*nb_obj:3]
zi = init[2:3*nb_obj:3]

# Figure initialization
fig = plt.figure(figsize=(figsize,figsize))
ax = fig.add_subplot(projection='3d',fc='black',box_aspect=(1,1,1))
sc = ax.scatter3D(xi,yi,zi,s=1.2,c='white')
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

# Animation
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

nframes = int(nsteps / sampling)
ani = FuncAnimation(fig, update,frames=nframes, interval=30, repeat=1) 

# Exporting Video
if save_vid:
    ani.save('anim.mp4',writer='ffmpeg',fps=18)
else:
    plt.show()