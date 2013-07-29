# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
import math
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
def contrastMultipleScatters(vec_r,particle_positions, wave_momentum):
        
        num_praticles = particle_positions[0]
                
        r = np.linalg.norm(vec_r)
        particle = 1
        
        A_holo = np.exp(-1j*wave_momentum*r)/r
        particle = 1
        while ( num_particles >= particle):
            r_ri = np.linalg.norm( vec_r - particle_positions[particle])
            A_holo += np.exp ( - 1j*wave_momentum * r_ri)/r_ri
            particle += 1
            
        #ans =  A_holo*np.conjugate(A_holo) - 1.0/(r*r) 
        ans = pow(np.abs(A_holo),2) - 1.0/(r*r)
 
        return ans

# <codecell>


# <codecell>

def contrastMatrix(particle_positions):
        vec_r = np.zeros(3)
        #vec_r1 = np.array( [0.,d_separation/2, z_obj])
        #vec_r2 = np.array( [0.,-d_separation/2, z_obj])
        #vec_r1 = np.array( vec1 )
        #vec_r2 = np.array( vec2 )
        #print "vec_r1: ", vec_r1
        #print "vec_r2: ", vec_r2
        I_mat = np.ndarray( (n_x_pixel, n_y_pixel))

        vec_r[2] = d_screen

        for i in range(n_x_pixel):
                vec_r[0] = x_arr[i]
                for j in range(n_y_pixel):
                        vec_r[1] = y_arr[j]
                        #print i, j
                        I_mat[i,j] = contrastMultipleScatters( vec_r, particle_positions, wave_momentum)  

                        
        plt.imshow(I_mat, extent = [-screen_width/2.0,screen_width/2.0,-screen_width/2.0,screen_width/2.0])
        plt.title('Contrast Image')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.colorbar()
        plt.show() 
        return I_mat

# <codecell>

def Kr(d_mesh,  nx, ny, I_mat, x_arr, y_arr, wave_momentum, d_screen):

#       d_mesh = -1 * d_mesh;
        global mesh_x_pix, mesh_y_pix, mesh_width
        si_norm_mat = np.zeros( (len(x_arr), len(y_arr), 3) )
        si_norm_mat[:,:,2] = d_screen

        for i in range(len(x_arr)):
                si_norm_mat[i,:,0] = x_arr[i]
        for j in range(len(y_arr)):
                si_norm_mat[:,j,1] = y_arr[j]

        vec_r = np.zeros(3)
#       x_mesh = np.arange(-nx*dx_mesh/2, nx*dx_mesh/2, dx_mesh)
#       y_mesh = np.arange(-ny*dy_mesh/2, ny*dy_mesh/2, dy_mesh)

        x_mesh = np.arange(nx)*mesh_x_pix - mesh_width/2.0
        y_mesh = np.arange(ny)*mesh_y_pix - mesh_width/2.0

        Kr_mesh = np.zeros( (len(x_mesh), len(y_mesh)), dtype = complex)

#       print "new mesh array"
        print "Kr func: ", len(x_mesh)
        for i in range(len(x_arr)):
                for j in range(len(y_arr)):
                        si_norm_mat[i,j,:] /= np.linalg.norm(si_norm_mat[i,j,:])

        vec_r[2] = d_mesh

        for i  in range(len(x_mesh)):
                vec_r[0] = x_mesh[i]
                for j in range(len(y_mesh)):
                        vec_r[1] = y_mesh[j]
                        for l in range(len(x_arr)):
                                for m in range(len(y_arr)):
                                        #theta = wave_momentum* np.dot(vec_r, si_norm_mat[l,m])
                                        #Kr_mesh[i,j] += I_mat[l,m]* np.exp(1j*theta)
                                        
                                        Kr_mesh[i,j] += I_mat[l,m]* np.exp(1j*wave_momentum*np.dot(vec_r, si_norm_mat[l,m]))
                                        
                                        #Kr_mesh[i,j] += I_mat[l,m]*math.cos(wave_momentum* np.dot(vec_r, si_norm_mat[l,m]))
                                        #Kr_mesh[i,j] += I_mat[l,m]*math.sin(wave_momentum* np.dot(vec_r, si_norm_mat[l,m]))*1j

                                        #Kr1_mesh[i,j] = abs(Kr_mesh[i,j])
                print i

        return Kr_mesh, x_mesh, y_mesh

# <codecell>

#pixel_pitch = 5
screen_width = 5000.0
#d_screen = screen_width/2 * math.sqrt(3.)
d_screen = 3500.0
print "Numerical Apperture is: ", (screen_width/2.)/ np.sqrt( pow((screen_width/2.),2) + pow(d_screen,2) )


z_obj = 120.0
n_x_pixel = 80
n_y_pixel = 80
size_x_pixel = screen_width/n_x_pixel
size_y_pixel = screen_width/n_y_pixel

d_separation = 0.8
wave_length = 0.4730
wave_momentum = 2*math.pi/wave_length

x_arr = np.arange(n_x_pixel)*size_x_pixel - screen_width/2.0
y_arr = np.arange(n_y_pixel)*size_y_pixel - screen_width/2.0

# <codecell>

num_particles = 2
particle_positions = [num_particles,   np.array( [ 0, d_separation/2, z_obj]), np.array([  1,-d_separation/2, z_obj])]
print particle_positions

# <codecell>

I_tilda = contrastMatrix( particle_positions )
#print I_tilda
print particle_positions

# <codecell>

contrastMultipleScatters(np.array([ 0., 0., d_screen]), particle_positions, wave_momentum)

# <codecell>

#for the object reconstruction
nx = 20
ny = 20
mesh_width = 10.0
mesh_x_pix = mesh_width / nx            #x-size of pixel of reconstruction mesh
mesh_y_pix = mesh_width / ny             #y-size of pixel of reconstruction mesh 
kr_120,x_mesh,y_mesh = Kr(z_obj,  nx, ny, I_tilda, x_arr, y_arr, wave_momentum, d_screen)

plt.imshow(np.abs(kr_120) , extent = [-mesh_width/2.0,mesh_width/2.0, -mesh_width/2.0 , mesh_width/2.0])
plt.title('Reconstruction')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.show()

# <codecell>

maxval = np.max(np.abs(kr_120))
plt.imshow(np.abs(kr_120)/maxval , extent = [-mesh_width/2.0,mesh_width/2.0, -mesh_width/2.0 , mesh_width/2.0])
plt.title('Reconstruction')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.show()

