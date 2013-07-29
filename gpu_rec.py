#librraries
import math
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import cv2.cv as cv
import cv2,sys
import time,math
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import scipy as sp
from numpy.linalg import *
from matplotlib.colors import LightSource
from matplotlib.backends.backend_pdf import PdfPages
import pycuda.cumath


screen_width=5000
screen_height=5000
z_screen=3500                   #distance between light source and screen
n_x_pixel=320                   #multiple of 80
n_y_pixel=320                   #multiple of 80
size_x_pixel=screen_width/n_x_pixel
size_y_pixel=screen_width/n_y_pixel

z_obj=50			#distance of particle from light source less then equal to 80
d_sep=0.8
wavelength=0.4730
wave_momentum=2*math.pi/wavelength
x_arr = np.arange(n_x_pixel)*size_x_pixel - screen_width/2.0
y_arr = np.arange(n_y_pixel)*size_y_pixel - screen_width/2.0

mesh_x_pixel=20
mesh_y_pixel=20
mesh_width=10.0
mesh_height=10.0
mesh_size_x=mesh_width/mesh_x_pixel
mesh_size_y=mesh_width/mesh_y_pixel
x_mesh=np.arange(mesh_x_pixel)*mesh_size_x - mesh_width/2
y_mesh=np.arange(mesh_y_pixel)*mesh_size_y - mesh_height/2
particlelist=[]
class particles:
	xpos=0.0;
	ypos=0.0;
	zpos=0.0;
	def __init__(self,x,y,z):
		self.xpos=x;
		self.ypos=y;
		self.zpos=z;

#specifying position of particles
#p1=particles(-1.0,2.0,z_obj)
#p2=particles(-1.0,-1.4,z_obj)
#p3=particles(2.0,-2.0,z_obj)
#p4=particles(-2.0,1.4,z_obj)
#p5=particles(1.4,0.4,z_obj)
#p6=particles(-1.0,-0.4,z_obj)
#p7=particles(-3.0,-3.0,z_obj)
#particlelist.append(p1)
#particlelist.append(p2)
#particlelist.append(p3)
#particlelist.append(p4)
#particlelist.append(p5)
#particlelist.append(p6)
radius=1.5
div=10000  #int(1000*1.2/radius)
theta=2*math.pi/div
for i in range(0,div):
	p=particles(radius*np.cos(i*theta),radius*np.sin(i*theta),z_obj)
	particlelist.append(p)
#particlelist.append(p7)
holograph = np.ndarray( (n_x_pixel , n_y_pixel ))
holograph = holograph.astype(np.complex64)
holograph2 = np.ndarray( (n_x_pixel , n_y_pixel ))
holograph2 = holograph2.astype(np.float32)
Kr_mesh = np.zeros( (len(x_mesh), len(y_mesh)), dtype = complex64)
j=np.complex64(-1j)
j=j*wave_momentum
Kr_real = np.ndarray((len(x_mesh),len(y_mesh)))
Kr_real = Kr_real.astype(np.float32)
Kr_img = np.ndarray((len(x_mesh),len(y_mesh)))
Kr_img = Kr_img.astype(np.float32)
def construct_holo():
	vec_R=np.zeros(3)
	vec_R[2]=z_screen

	for i in range(n_x_pixel):
		vec_R[0]=x_arr[i]
		for j in range(n_y_pixel):
			vec_R[1]=y_arr[j]
			holograph[i,j] = 0.0
			holograph[i,j] = constructmultipleparticles( vec_R)
	plot_holo(holograph)
def plot_holo(holograph):
	holograph1 = np.ndarray( (n_x_pixel , n_y_pixel ))
	vec_R=np.zeros(3)
	vec_R[2]=z_screen
	for i in range(n_x_pixel):
		vec_R[0]=x_arr[i]
		for j in range(n_y_pixel):
			vec_R[1]=y_arr[j]
			mod_R=np.linalg.norm(vec_R)
			holograph[i,j]=holograph[i,j]+np.exp(-1j*wave_momentum*mod_R)/mod_R
			holograph1[i,j]=pow(np.abs(holograph[i,j]),2) - 1.0/(mod_R*mod_R)
			
	plt.imshow(holograph1, extent = [-screen_width/2.0,screen_width/2.0,-screen_width/2.0,screen_width/2.0])
	plt.title('Contrast Image')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.colorbar()
	plt.show()
	return holograph1
def constructmultipleparticles(vec_R):
	mod_R=np.linalg.norm(vec_R)
	A_holo=0
#	A_holo = np.exp(-1j*wave_momentum*mod_R)/mod_R
	particles=0
	no_particles=len(particlelist);
	vec_r=np.zeros(3)
	vec_r[2]=z_obj
	while(no_particles>particles):
		vec_r[0]=particlelist[particles].xpos
		vec_r[1]=particlelist[particles].ypos
		r_ri = np.linalg.norm( vec_R - vec_r)
		A_holo +=np.exp(j * r_ri)/r_ri
		particles =particles + 1
	return A_holo
#	ans = pow(np.abs(A_holo),2) - 1.0/(mod_R*mod_R)
#	return ans

def reconstruction(holograph):
	vec_r=np.zeros(3)
	vec_R=np.zeros(3)
	vec_r[2]=z_obj
	vec_R[2]=z_screen
#	Kr_mesh = np.zeros( (len(x_mesh), len(y_mesh)), dtype = complex)
	for i in range(len(x_mesh)):
		vec_r[0]=x_mesh[i]
		for j in range(len(y_mesh)):	
			vec_r[1]=y_mesh[j]
			for l in range(len(x_arr)):
				vec_R[0]=x_arr[l]
				for m in range(len(y_arr)):
					vec_R[1]=y_arr[m]
					R_hat=np.linalg.norm(vec_R)
					Kr_mesh[i,j]+= holograph[l,m]* np.exp(1j*wave_momentum*(np.dot(vec_r, vec_R)/R_hat))
	print Kr_mesh.real
	plt.imshow(np.abs(Kr_mesh) , extent = [-mesh_width/2.0,mesh_width/2.0, -mesh_width/2.0 , mesh_width/2.0])
	plt.title('Reconstruction')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.colorbar()
	plt.show()

#construct_holo()
#reconstruction()


#parallel code
block_x=mesh_x_pixel
block_y=mesh_y_pixel
block_1=1
grid_x=int(math.ceil(float(n_x_pixel)/block_x))
grid_y=int(math.ceil(float(n_y_pixel)/block_y))
grid_z=int(len(particlelist))

#convering to gpu format
holograph = holograph.astype(np.complex64)
Kr_mesh=Kr_mesh.astype(np.complex64)
Kr_mesh_gpu=cuda.mem_alloc(Kr_mesh.nbytes)
Kr_real_gpu=cuda.mem_alloc(Kr_real.nbytes)
Kr_img_gpu=cuda.mem_alloc(Kr_img.nbytes)
hg_gpu = cuda.mem_alloc(holograph.nbytes)
cuda.memcpy_htod(hg_gpu, holograph)
cuda.memcpy_htod(Kr_mesh_gpu,Kr_mesh)
cuda.memcpy_htod(Kr_real_gpu,Kr_real)
cuda.memcpy_htod(Kr_img_gpu,Kr_img)
particle_array = np.ndarray((len(particlelist),3))
particle_array = particle_array.astype(np.float32)
for i in range(len(particlelist)):
	particle_array[i,0]=particlelist[i].xpos
	particle_array[i,1]=particlelist[i].ypos
	particle_array[i,2]=particlelist[i].zpos
#particle_array = particle_array.astype(np.float32)
particle_array_gpu=cuda.mem_alloc(particle_array.nbytes)
cuda.memcpy_htod(particle_array_gpu,particle_array)
z_screen_gpu = np.int32(z_screen)
wave_momentum_gpu = np.float32(wave_momentum)
screen_width_gpu=np.float32(screen_width)
screen_height_gpu=np.float32(screen_height)
size_x_pixel_gpu=np.float32(size_x_pixel)
size_y_pixel_gpu=np.float32(size_y_pixel)
mesh_width_gpu=np.float32(mesh_width)
mesh_height_gpu=np.float32(mesh_height)
mesh_size_x_pixel_gpu=np.float32(mesh_size_x)
mesh_size_y_pixel_gpu=np.float32(mesh_size_y)
j_gpu=np.complex64(j)
z_obj_gpu=np.int32(z_obj)
ans=1.0
ans_gpu=np.float32(ans)
#include <pycuda-complex.hpp>
mod = SourceModule("""
#include <pycuda-complex.hpp>
__global__ void constructHolo(pycuda::complex<float> *hg,float *particleArray, int zScreen,
                        float sizeXPixel,float sizeYPixel,float screenWidth , float screenHeight ,  float wavemomentum , pycuda::complex<float> j_gpu)
  {
    int idx = (blockIdx.x* blockDim.x + threadIdx.x);            //iteration over x direction of grid
    int idy = (blockIdx.y* blockDim.y + threadIdx.y);            //iteration over y direction of grid
    int idz = (blockIdx.z* blockDim.z + threadIdx.z);            //iteration over z direction of grid
    int index = idx + idy*blockDim.x*gridDim.x;                  //index generation of 2D matrix in grid
    
    if(idx < screenWidth/sizeXPixel ){
    
        //vector to a point on holograph (R)
        //fmod required as idx may become greater than nXPixels(number of pixels in X direction)
        float xR = fmod((float)idx,(float)(screenWidth/sizeXPixel)) *sizeXPixel - screenWidth/2.0; 
        float yR = fmod((float)idy,(float)(screenHeight/sizeYPixel))*sizeYPixel - screenHeight/2.0;
        float zR = zScreen; 
        float normR = sqrt(pow(xR,2) + pow(yR,2) + pow(zR,2));
        //screenWidth/sizeXPixel gives number of pixels on screen in x direction (nXPixels)
    
        //vector to a particle (r)
        float xr = particleArray[(idz*3)];                                 
        float yr = particleArray[(idz*3) + 1];
        float zr = particleArray[(idz*3) + 2];
        float normr = sqrt(pow(xr,2) + pow(yr,2) + pow(zr,2));
        float normRr = sqrt(pow(xR-xr,2) + pow(yR-yr,2) + pow(zR-zr,2));            //norm of (R - r) 
	hg[index] +=(exp(j_gpu * (normRr) ) ) / (normRr );

   }
   }
		""")
mod2 = SourceModule("""
#include <pycuda-complex.hpp>
__global__ void reconstruct_holo(float *Kr_real,float *Kr_img,float *hg , int z_screen  , int z_obj , float wavemomentum , pycuda::complex<float> j_gpu , float screen_width , float screen_height , float size_x_pixel , float size_y_pixel , float mesh_height , float mesh_width , float mesh_size_x_pixel ,
	float mesh_size_y_pixel)
{
	int index_mesh=threadIdx.x + threadIdx.y*blockDim.x;
	int index_holo=blockIdx.x + blockIdx.y*gridDim.x;
	int k=threadIdx.x;
	int l=threadIdx.y;
	int i=blockIdx.x;
	int j=blockIdx.y;
//	if(index_holo<6400 && index_mesh<400 )

//	{
		//vector to point on holograph
		float xR=fmod((float)i,(float)(screen_width/size_x_pixel)) *size_x_pixel - screen_width/2.0; 
		float yR=fmod((float)j,(float)(screen_height/size_y_pixel))*size_y_pixel - screen_height/2.0;
		float zR=z_screen;
		float xr=fmod((float)k,(float)(mesh_width/mesh_size_x_pixel)) *mesh_size_x_pixel - mesh_width/2.0; 
		float yr=fmod((float)l,(float)(mesh_height/mesh_size_y_pixel))*mesh_size_y_pixel - mesh_height/2.0;
		float zr=z_obj;
	        float normR = sqrt(pow(xR,2) + pow(yR,2) + pow(zR,2));
		int dot = xr*xR + yr*yR + zr*zR;
		atomicAdd(&Kr_real[index_mesh],(hg[index_holo]*(cos((wavemomentum*dot/normR)))));
		atomicAdd(&Kr_img[index_mesh],(hg[index_holo]*(sin((wavemomentum*dot/normR)))));
//	}


}

""")

def construct_holo_gpu():
#    j=np.complex64(-1j)
#    j=j*wave_momentum
#    print j
#    j_gpu=np.complex64(j)
    func = mod.get_function("constructHolo")
    kernel = func(hg_gpu,particle_array_gpu, z_screen_gpu,size_x_pixel_gpu,size_y_pixel_gpu,screen_width_gpu,screen_height_gpu,
                  wave_momentum_gpu,j_gpu, 
                  block=(block_x,block_y,1),
                  grid=(grid_x,grid_y,grid_z),
                  time_kernel = True)
#    print holograph
    cuda.memcpy_dtoh(holograph,hg_gpu)
#    print holograph
    holograph2=plot_holo(holograph)
#    plt.imshow(holograph2, extent = [-screen_width/2.0,screen_width/2.0,-screen_width/2.0,screen_width/2.0])
#    plt.title('Contrast Image')
#    plt.xlabel('x')
#    plt.ylabel('y')
#    plt.colorbar()
#    plt.show()
    holograph2 = holograph2.astype(np.float32)
    hg2_gpu = cuda.mem_alloc(holograph2.nbytes)
    cuda.memcpy_htod(hg2_gpu, holograph2)
    func=mod2.get_function("reconstruct_holo")
    kernel = func(Kr_real_gpu,Kr_img_gpu,hg2_gpu,z_screen_gpu,z_obj_gpu,wave_momentum_gpu,-j_gpu,screen_width_gpu,screen_height_gpu,
		    size_x_pixel_gpu,size_y_pixel_gpu,mesh_height_gpu,mesh_width_gpu,mesh_size_x_pixel_gpu,
		    mesh_size_y_pixel_gpu,
		    block=(mesh_x_pixel,mesh_y_pixel,1),
		    grid=(n_x_pixel,n_y_pixel,1),
		    time_kernel = True)
    cuda.memcpy_dtoh(Kr_real,Kr_real_gpu)
    cuda.memcpy_dtoh(Kr_img,Kr_img_gpu)
    cuda.memcpy_dtoh(holograph2,hg2_gpu)
#    print Kr_real
    Kr_mesh_gpu.free()
    max1=0.0
    for i in range(mesh_x_pixel):
	    for j in range(mesh_y_pixel):
		    Kr_mesh[i,j]=Kr_real[i,j] + 1j*Kr_img[i,j]
		    Kr_real[i,j]=np.abs(Kr_mesh[i,j])
		    if(max1<Kr_real[i,j]):
			    max1=Kr_real[i,j]
    for i in range(mesh_x_pixel):
	    for j in range(mesh_y_pixel):
		    if (Kr_real[i,j]< 14*max1/16):
			    Kr_real[i,j]=0.0
    plt.imshow(Kr_real , extent = [-mesh_width/2.0,mesh_width/2.0, mesh_width/2.0 , -mesh_width/2.0])
    plt.title('Reconstruction')
    plt.xlabel('x')
    plt.gca().invert_yaxis()
    plt.ylabel('y')
    plt.colorbar()
    plt.show()

#    reconstruction(holograph2)

construct_holo_gpu()
