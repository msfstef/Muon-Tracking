# S. Zenz, Imperial College
# April 2015
# A starting script to draw the CMS solenoid in 2d and 3d and plot a track and hits

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# Set to false to skip making 3D plot
do3D = False

# Features of CMS detector
B = 3.8 # Magnetic field strength, units: Tesla
r_sol = 2.950 #~ Radius of solenoid volume, units: meters
r_muon = r_sol + 3. # Radius of muon ironman tracker, units: meters
z_sol = 3.00 # Length of solenoid volume, units: meters

# Tracker layers (z = 0) [units: m]
# These are approximate and read from http://www.hephy.at/user/friedl/diss/html/node26.html
# You can see there that the layout is rather complicated, but let's keep things simple for now
r_layers = [0.04,0.07,0.11,0.26,0.32,0.43,0.62,0.71,0.79,0.88,0.97,1.07]

# initial particle momentum
pz = 0. # units: GeV/c
pt = 3. # units: GeV/c

# initial particle position
x0 = 0. # units: meters
y0 = 0. # units: meters
z0 = 0. # units: meters


# initial angle in x-y plane
a = 0
phi0 = a*(np.pi/2) - np.pi/2 # 90 degree conventional correction

# Sense of direction determined by charge (q = -1 muon, q = 1 antimuon)
q = -1

# Create 3D axes
if do3D:
  fig3 = plt.figure(1)
  ax3 = fig3.add_subplot(111, projection='3d')
  ax3.set_xlabel("X [m]")
  ax3.set_ylabel("Y [m]")
  ax3.set_zlabel("Z [m]")
  ax3.set_xlim3d(-1.1*r_muon,1.1*r_muon)
  ax3.set_ylim3d(-1.1*r_muon,1.1*r_muon)
  ax3.set_zlim3d(-1.1*r_muon,1.1*r_muon)

# Create 2d axes
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.set_xlabel("X [m]")
ax2.set_ylabel("Y [m]")
ax2.set_xlim(-1.1*r_muon,1.1*r_muon)
ax2.set_ylim(-1.1*r_muon,1.1*r_muon)

# Create solenoid points
x_grid=np.linspace(-r_sol, r_sol, 100)
z_grid=np.linspace(-z_sol, z_sol, 100)
Xc, Zc=np.meshgrid(x_grid, z_grid)
Yc = np.sqrt(r_sol*r_sol-Xc*Xc)
y_grid = np.sqrt(r_sol*r_sol-x_grid*x_grid)

#~ Create 2nd solenoid points
x_grid2=np.linspace(-r_muon, r_muon, 200)
z_grid2=np.linspace(-z_sol, z_sol, 200)
Xc2, Zc2=np.meshgrid(x_grid2, z_grid2)
Yc2 = np.sqrt(r_muon*r_muon-Xc2*Xc2)
y_grid2 = np.sqrt(r_muon*r_muon-x_grid2*x_grid2)


# Plot solenoid points in 2d
ax2.plot(x_grid,y_grid,color='b',label='CMS solenoid')
ax2.plot(x_grid,-y_grid,color='b')

#~ Plot solenoid 2 points in 2d
ax2.plot(x_grid2,y_grid2,color='c',label='Muon solenoid')
ax2.plot(x_grid2,-y_grid2,color='c')


# Plot solenoid points as mesh in 3d
if do3D:
  ax3.plot_surface(Xc, Yc, Zc,  rstride=4, cstride=4, color='b', alpha=0.2)
  ax3.plot_surface(Xc, -Yc, Zc,  rstride=4, cstride=4, color='b', alpha=0.2)
  
# Plot solenoid 2 points as mesh in 3d
if do3D:
  ax3.plot_surface(Xc2, Yc2, Zc2,  rstride=4, cstride=4, color='c', alpha=0.2)
  ax3.plot_surface(Xc2, -Yc2, Zc2,  rstride=4, cstride=4, color='c', alpha=0.2)

# Derivations of motion in magnetic field and helix parameterization
# http://www.physics.iitm.ac.in/~sercehep2013/track2_Gagan_Mohanty.pdf

# Calculate radius of curvature
Rc = pt / (0.3*B) # units: meters

# Calculate dip angle
ldip = np.arctan(pz/pt)

# Create parameters of helix
s = np.linspace(0, 16*np.pi*Rc, 100000)  # path length
z = s*np.sin(ldip)  + z0
x = x0 -q*Rc*(np.cos(phi0 -q*s*np.cos(ldip)/Rc) - np.cos(phi0))
y = y0 -q*Rc*(np.sin(phi0 -q*s*np.cos(ldip)/Rc) - np.sin(phi0))



xhit = []
yhit = []
zhit = []

# don't plot helix beyond tracker volume
r = np.sqrt(x*x+y*y)
for i in range(len(r)):
  for rl in r_layers:
    if i>0 and r[i] > rl and r[i-1] < rl:
      xhit.append(0.5*(x[i]+x[i-1]))
      yhit.append(0.5*(y[i]+y[i-1]))
      zhit.append(0.5*(z[i]+z[i-1]))
  if r[i] > r_sol or abs(z[i]) > z_sol:
    print "Crossing tracker at %ith point which has (r,z)=(%f,%f) (x,y)=(%f,%f) z/r=%f pz/pt=%f Rc=%f" % (i,r[i],z[i],x[i],y[i],z[i]/r[i],pz/pt,Rc)
    crosspoint = i
    x_fin = x[i]
    y_fin = y[i]
    z_fin = z[i]
    x=x[:i]
    y=y[:i]
    z=z[:i]
    break
    
    

    
#~ Create parameters of second part of helix
s = np.linspace(0, 16*np.pi*Rc, 100000)  # path length
tangent_up = (x0-x_fin + q*Rc*np.cos(phi0))
tangent_down = (y_fin - y0 - q*Rc*np.sin(phi0))
if q > 0 :
    phi1 = np.arctan2(tangent_up,tangent_down) + 3*np.pi/2
else:
    phi1 = np.arctan2(tangent_up,tangent_down) + np.pi/2
    
z2 = s*np.sin(ldip)  + z_fin
x2 = x_fin + q*Rc*(np.cos(phi1 +q*s*np.cos(ldip)/Rc) - np.cos(phi1))
y2 = y_fin + q*Rc*(np.sin(phi1 +q*s*np.cos(ldip)/Rc) - np.sin(phi1))

#~ don't plot second helix beyond muon solenoid volume
r2 = np.sqrt(x2*x2+y2*y2)
for i in range(len(r2)):
  #~ for rl in r_layers:
    #~ if i>0 and r[i] > rl and r[i-1] < rl:
      #~ xhit.append(0.5*(x[i]+x[i-1]))
      #~ yhit.append(0.5*(y[i]+y[i-1]))
      #~ zhit.append(0.5*(z[i]+z[i-1]))
  if r2[i] > r_muon or abs(z2[i]) > z_sol:
    print "Truncating at %ith point which has (r,z)=(%f,%f) (x,y)=(%f,%f) z/r=%f pz/pt=%f Rc=%f" % (i+crosspoint,r2[i],z2[i],x2[i],y2[i],z2[i]/r2[i],pz/pt,Rc)
    x2=x2[:i]
    y2=y2[:i]
    z2=z2[:i]
    break




# Plot the helix and hits
if do3D: ax3.plot(np.append(x,x2),np.append(y,y2),np.append(z,z2), label='Particle path', color='r')
if do3D: ax3.plot(xhit,yhit,zhit, label="Tracker hit (unsmeared)",color='r',marker='o',linestyle="none")
ax2.plot(np.append(x,x2),np.append(y,y2),label='Particle path',color='r')
ax2.plot(np.array(xhit),np.array(yhit),label="Tracker hit (unsmeared)",color='r',marker='o',linestyle="none")
#~ ax2.legend(numpoints=1)

# Show the figures we've made
plt.show()
