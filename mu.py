# Initial Version Created by S. Zenz, Imperial College
# Developed by S. Mousafeiris and N. Koukoulekidis, Imperial College
# May 2015
# A starting script to draw the CMS solenoid in 2d and 3d and plot a track and hits


import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# Features of CMS detector
B = 3.8 # Magnetic field strength in 1st solenoid, units: Tesla
B2 = 2.0 # Magnetic field strength in 2nd solenoid, units: Tesla
r_sol = 2.950 #~ Radius of solenoid volume, units: meters
r_sol2 = r_sol + 3.550 # Radius of muon ironman tracker, units: meters
z_sol = 3.00 # Length of solenoid volume, units: meters

# Tracker layers (z = 0) [units: m]
# These are approximate and read from http://www.hephy.at/user/friedl/diss/html/node26.html
# You can see there that the layout is rather complicated, but let's keep things simple for now
r_layers = [0.04,0.07,0.11,0.26,0.32,0.43,0.62,0.71,0.79,0.88,0.97,1.07]

#~ # initial particle momentum
#~ pz = momentum in z direction # units: GeV/c
#~ pt  = transverse momentum # units: GeV/c

#~ # initial particle position
#~ x0 = initial x position # units: meters
#~ y0 = initial y position # units: meters
#~ z0 = initial z position # units: meters


#~ # initial angle in x-y plane
#~ phi0 = initial angle of particle in transverse plane # units : rads

#~ # Sense of direction determined by charge 
#~ q = -1 muon (anticlockwise) , q = 1 antimuon (clockwise)

#Plots particle path and tracker hits in crosssection of detector.
def plot_tracker_path(x_path, y_path, z_path, xhit, yhit, zhit):
	# Create 2d axes
	fig2 = plt.figure(2)
	ax2 = fig2.add_subplot(111)
	ax2.set_xlabel("X [m]")
	ax2.set_ylabel("Y [m]")
	ax2.set_xlim(-1.1*r_sol2,1.1*r_sol2)
	ax2.set_ylim(-1.1*r_sol2,1.1*r_sol2)
	
	solenoid_points(ax2)
	plot_path(x_path, y_path, z_path, xhit, yhit, zhit, ax2)
	

#Plots particle path and tracker hits in detector, 3D plot.	
def plot_tracker_path_3D(x_path, y_path, z_path, xhit, yhit, zhit):
	# Create 3D axes
	fig3 = plt.figure(1)
	ax3 = fig3.add_subplot(111, projection='3d')
	ax3.set_xlabel("X [m]")
	ax3.set_ylabel("Y [m]")
	ax3.set_zlabel("Z [m]")
	ax3.set_xlim3d(-1.1*r_sol2,1.1*r_sol2)
	ax3.set_ylim3d(-1.1*r_sol2,1.1*r_sol2)
	ax3.set_zlim3d(-1.1*z_sol,1.1*z_sol)

	solenoid_points_3D(ax3)
	plot_path_3D(x_path, y_path, z_path, xhit, yhit, zhit, ax3)

# Creates the detector in 2D.
def solenoid_points(ax2):
	# Create solenoid points
	x_grid=np.linspace(-r_sol, r_sol, 100)
	y_grid = np.sqrt(r_sol*r_sol-x_grid*x_grid)

	# Create 2nd solenoid points
	x_grid2=np.linspace(-r_sol2, r_sol2, 200)
	y_grid2 = np.sqrt(r_sol2*r_sol2-x_grid2*x_grid2)

	# Plot solenoid points in 2d
	ax2.plot(x_grid,y_grid,color='b',label='CMS solenoid')
	ax2.plot(x_grid,-y_grid,color='b')

	# Plot 2nd solenoid points in 2d
	ax2.plot(x_grid2,y_grid2,color='c',label='Muon solenoid')
	ax2.plot(x_grid2,-y_grid2,color='c')

#Creates the detector in 3D.		
def solenoid_points_3D(ax3):
	# Create solenoid points
	x_grid=np.linspace(-r_sol, r_sol, 100)
	z_grid=np.linspace(-z_sol, z_sol, 100)
	Xc, Zc=np.meshgrid(x_grid, z_grid)
	Yc = np.sqrt(r_sol*r_sol-Xc*Xc)
	y_grid = np.sqrt(r_sol*r_sol-x_grid*x_grid)

	# Create 2nd solenoid points
	x_grid2=np.linspace(-r_sol2, r_sol2, 200)
	z_grid2=np.linspace(-z_sol, z_sol, 200)
	Xc2, Zc2=np.meshgrid(x_grid2, z_grid2)
	Yc2 = np.sqrt(r_sol2*r_sol2-Xc2*Xc2)
	y_grid2 = np.sqrt(r_sol2*r_sol2-x_grid2*x_grid2)

	# Plot solenoid points as mesh in 3d
	ax3.plot_surface(Xc, Yc, Zc,  rstride=4, cstride=4, color='b', alpha=0.2)
	ax3.plot_surface(Xc, -Yc, Zc,  rstride=4, cstride=4, color='b', alpha=0.2)
	
	# Plot 2nd solenoid points as mesh in 3d
	ax3.plot_surface(Xc2, Yc2, Zc2,  rstride=4, cstride=4, color='c', alpha=0.2)
	ax3.plot_surface(Xc2, -Yc2, Zc2,  rstride=4, cstride=4, color='c', alpha=0.2)


#Calculates the particle's angle with the x-axis at the point where it crosses the first part of the detector.
def phi_calc(Rc, x_fin, y_fin, phi0, x0, y0, q):
    	tangent_up = (x0-x_fin + q*Rc*np.cos(phi0))
	tangent_down = (y_fin - y0 - q*Rc*np.sin(phi0))
	if q > 0 :
		phi1 = np.arctan2(tangent_up,tangent_down) + 3*np.pi/2
	else:
		phi1 = np.arctan2(tangent_up,tangent_down) + np.pi/2
	return phi1


#Calculates the path of the particle through a layer of the detector and the tracker hits.
def path_creation(Rc, ldip, x0, y0, z0, phi0, q, pt, pz, layer, xhit, yhit, zhit, x_path, y_path, z_path, crosspoint):
	# Create parameters of helix
	s = np.linspace(0, 8*np.pi*Rc, 100000)  # path length
	z = s*np.sin(ldip)  + z0
	x = x0 -q*Rc*(np.cos(phi0 -q*s*np.cos(ldip)/Rc) - np.cos(phi0))
	y = y0 -q*Rc*(np.sin(phi0 -q*s*np.cos(ldip)/Rc) - np.sin(phi0))
	# don't plot helix beyond tracker volume
	r = np.sqrt(x*x+y*y)
	for i in range(len(r)):
		for rl in r_layers:
			if i>0 and r[i] > rl and r[i-1] < rl:
				xhit.append(0.5*(x[i]+x[i-1]))
				yhit.append(0.5*(y[i]+y[i-1]))
				zhit.append(0.5*(z[i]+z[i-1]))
		
		if layer == 0 and r[i] > r_sol:
			print "Crossing tracker at %ith point which has (r,z)=(%f,%f) (x,y)=(%f,%f) z/r=%f pz/pt=%f Rc=%f" % (i,r[i],z[i],x[i],y[i],z[i]/r[i],pz/pt,Rc)
			crosspoint = i
			x_fin = x[i]
			y_fin = y[i]
			z_fin = z[i]
			x_path= np.append(x_path, x[:i])
			y_path= np.append(y_path, y[:i])
			z_path= np.append(z_path, z[:i])
			return x_fin, y_fin, z_fin, x_path, y_path, z_path, crosspoint
			break
		if layer == 0 and abs(z[i]) > z_sol:
		    print "Truncating at %ith point which has (r,z)=(%f,%f) (x,y)=(%f,%f) z/r=%f pz/pt=%f Rc=%f" % (i,r[i],z[i],x[i],y[i],z[i]/r[i],pz/pt,Rc)
                    crosspoint = i
                    x_fin = x[i]
                    y_fin = y[i]
                    z_fin = z[i]
                    x_path= np.append(x_path, x[:i])
                    y_path= np.append(y_path, y[:i])
                    z_path= np.append(z_path, z[:i])
                    return x_fin, y_fin, z_fin, x_path, y_path, z_path, crosspoint
                    break
		
		if layer == 1 and  r[i] > r_sol2 or abs(z[i]) > z_sol:
		    print "Truncating at %ith point which has (r,z)=(%f,%f) (x,y)=(%f,%f) z/r=%f pz/pt=%f Rc=%f" % (i+crosspoint,r[i],z[i],x[i],y[i],z[i]/r[i],pz/pt,Rc)
		    x_path= np.append(x_path, x[:i])
		    y_path= np.append(y_path, y[:i])
		    z_path= np.append(z_path, z[:i])
		    return x_path, y_path, z_path, xhit, yhit, zhit
		    break
			
	
	
#Creates full, continuous path of particle throughout detector.
def particle_path(x0, y0, z0, phi0, q, pt, pz, xhit, yhit, zhit, x_path, y_path, z_path):
	# Calculate radius of curvature
	Rc = pt / (0.3*B) # units: meters
	Rc2 = pt / (0.3*B2) # units: meters
	# Calculate dip angle
	ldip = np.arctan(pz/pt)
	layer = 0
	x_fin, y_fin, z_fin, x_path, y_path, z_path, crosspoint = path_creation(Rc, ldip, x0, y0, z0, phi0, q, pt, pz, layer, xhit, yhit, zhit, x_path, y_path, z_path, crosspoint=0)
	if np.sqrt(x_fin*x_fin + y_fin*y_fin) < r_sol:
	    return x_path, y_path, z_path, xhit, yhit, zhit
	else :
	   layer = 1
	   phi1 = phi_calc(Rc, x_fin, y_fin, phi0, x0, y0, q)
	   x_path, y_path, z_path, xhit, yhit, zhit = path_creation(Rc2, ldip, x_fin, y_fin, z_fin, phi1, -q, pt, pz, layer, xhit, yhit, zhit, x_path, y_path, z_path, crosspoint)
	   return x_path, y_path, z_path, xhit, yhit, zhit

#Plots particle path and hits in cross-section of detector.
def plot_path(x_path, y_path, z_path, xhit, yhit, zhit, ax2):
	# Plot the helix and hits
	ax2.plot(x_path, y_path,label='Particle path',color='r')
	ax2.plot(np.array(xhit),np.array(yhit),label="Tracker hit (unsmeared)",color='r',marker='o',linestyle="none")
	ax2.legend(numpoints=1)

	# Show the figures we've made
	plt.show()
	
#Plots particle path and hits in 3D detector.
def plot_path_3D(x_path, y_path, z_path, xhit, yhit, zhit, ax3):
	# Plot the helix and hits
	ax3.plot(x_path, y_path, z_path, label='Particle path', color='r')
	ax3.plot(xhit,yhit,zhit, label="Tracker hit (unsmeared)",color='r',marker='o',linestyle="none")
	# Show the figures we've made
	plt.show()


#Combines all previous functions, takes initial conditions of particle as arguments and plots its path in desired way.
def muon_path(x0=0., y0=0., z0=0., phi0=0., pt=3., pz=1., q = -1, do3D=False):
	phi0 = phi0 -np.pi/2 # 90 degree conventional correction
	x_path,y_path,z_path,xhit,yhit,zhit = [], [], [], [], [], []
	x_path,y_path,z_path,xhit,yhit,zhit = particle_path(x0, y0, z0, phi0, q, pt, pz, xhit, yhit, zhit, x_path, y_path, z_path)
	plot_tracker_path(x_path, y_path, z_path, xhit, yhit, zhit)
	if do3D: plot_tracker_path_3D(x_path, y_path, z_path, xhit, yhit, zhit)
	
	
muon_path()
