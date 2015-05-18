from mu import *
from scipy import optimize

# Errors

err_1_t=10.*10**(-6) # Pixel detector transverse error, units: meters
err_1_z=10.*10**(-6) # Pixel detector z error, units: meters
err_1 = [err_1_t, err_1_z]

err_2_t=20.*10**(-6) # Internal Silicon Strip Tracker transverse error, units: meters
err_2_z=20.*10**(-6) # Internal Silicon Strip Tracker z error, units: meters
err_2 = [err_2_t, err_2_z]

err_3_t=30.*10**(-6) # External Silicon Strip Tracker transverse error, units: meters
err_3_z=30.*10**(-6) # External Silicon Strip Tracker z error, units: meters
err_3 = [err_3_t, err_3_z]

# Incorrect values used for testing
err_4_t=30.*10**(-6) # MB transverse error, units: meters
err_4_z=30.*10**(-6) # MB z error, units: meters
err_4 = [err_4_t, err_4_z]

r_1 = [0.1, 0.2, 0.3]
r_2 = [0.4, 0.5, 0.6]
r_3 = [0.7, 0.8, 0.9]
r_4 = [2. , 2.5, 2.8]


#Generates data (tracker hit locations) from particle with specified initial conditions.
def gen_data(x0=0., y0=0., z0=0., phi0=0., pt=3., pz=1., q=-1):
    xhit, yhit, zhit = np.asarray(muon_path(x0, y0, z0, phi0, pt, pz, q, plot = False))
    rhit = np.sqrt(xhit*xhit + yhit*yhit)
    phihit = np.arctan2(yhit, xhit)
    return rhit, phihit, zhit 


#Smears data by assigning each data point a value based on a normal distribution
#created using the tracking device's standard measuring error.
def smear_func(initial, final, err, rhit, phihit, zhit):
    i = initial
    f = i + len(final)
    r_smear = rhit[i:f]
    phi_smear = np.random.normal(phihit[i:f], err[0]/rhit[i:f])
    x_smear = r_smear * np.cos(phi_smear)
    y_smear = r_smear * np.sin(phi_smear)
    z_smear = np.random.normal(zhit[i:f], err[1])
    return f, x_smear, y_smear, z_smear

#Smears each data point in the appropriate way (data from multiple trackers).    
def smear_data(rhit, phihit, zhit):
    x_smear, y_smear, z_smear = [],[],[]
    err_layers = [err_1, err_2, err_3, err_4]
    det_layers = [[], r_1, r_2, r_3, r_4]
    temp = [0]
    for j in range(len(det_layers)-1):
        temp = smear_func(temp[0], det_layers[j+1], err_layers[j], rhit, phihit, zhit)
        x_smear = np.append(x_smear, temp[1])
        y_smear = np.append(y_smear, temp[2])
        z_smear = np.append(z_smear, temp[3])
    return x_smear, y_smear, z_smear


#The following functions are used to create a least squares circle fit.
#======================================================================#
def calc_R(x,y, xc, yc):
    # Calculate the distance of each 2D points from the centre (xc, yc).
    return np.sqrt((x-xc)**2 + (y-yc)**2)
 
def f(c, x, y):
    # Calculate the algebraic distance between the data points and the mean circle centred at c=(xc, yc).
    Ri = calc_R(x, y, *c)
    return Ri - Ri.mean()
 
def leastsq_circle(x,y):
    # Coordinates of the barycenter
    x_m = np.mean(x)
    y_m = np.mean(y)
    center_estimate = x_m, y_m
    center, ier = optimize.leastsq(f, center_estimate, args=(x,y))
    xc, yc = center
    Ri       = calc_R(x, y, *center)
    R        = Ri.mean()
    #R_err    = Ri.std()/np.sqrt(len(Ri))
    #residu   = np.sum((Ri - R)**2)
    return xc, yc, R
 
def plot_data_circle(x,y, xc, yc, R):
    f = plt.figure( facecolor='white')  #figsize=(7, 5.4), dpi=72,
    plt.axis('equal')
 
    theta_fit = np.linspace(-np.pi, np.pi, 180)
 
    x_fit = xc + R*np.cos(theta_fit)
    y_fit = yc + R*np.sin(theta_fit)
    plt.plot(x_fit, y_fit, 'b-' , label="fitted circle", lw=2)
    plt.plot([xc], [yc], 'bD', mec='y', mew=1)
    plt.xlabel('x')
    plt.ylabel('y')   
    # plot data
    plt.plot(x, y, 'r.', label='data', mew=1)
 
    plt.legend(loc='best',labelspacing=0.1 )
    plt.grid()
    plt.title('Least Squares Circle')
    plt.show()
#===================================================================#  

#Calculates transverse momentum by fitting circular path in smeared data.
def pt_calc(rhit, phihit, zhit):
    x_data, y_data, z_data = smear_data(rhit, phihit, zhit)
    xc, yc, Rc = leastsq_circle(x_data, y_data)
    pt = 0.3*B*Rc
    return pt

#Iterates calculation of transverse momentum over a given amount of times,
#and calculates a mean and a standard error for the estimated transverse
#momentum for this particle.
def pt_datapoint(pt, iter_num):
    rhit, phihit, zhit = gen_data(pt=pt)
    pt_data = []
    for i in range(iter_num):
        pt_data = np.append(pt_data, pt_calc(rhit, phihit, zhit))
    pt_mean = np.mean(pt_data)
    pt_err = (np.std(pt_data))/(np.sqrt(iter_num))
    return pt_mean, pt_err


#Iterates previous function over a given range of transverse momenta.
def pt_data(pt_i, pt_f, point_num, iter_num):
    pt, pt_actual, pt_err = [], [], []
    stepsize = float((pt_f - pt_i))/point_num
    for i in range(int((pt_f-pt_i)/stepsize)):
        temp1, temp2 = pt_datapoint(iter_num = iter_num, pt = pt_i + i*stepsize)
        pt = np.append(pt, temp1)
        pt_err = np.append(pt_err, temp2)
        pt_actual = np.append(pt_actual, pt_i +i*stepsize)
    return pt, pt_err, pt_actual
    
#Plots the standard error in the transverse momentum and its deviation from
#the actual value against the transverse momentum.
def plot_data(pt_i=1, pt_f=6, point_num=100, iter_num=200):
    pt, pt_err, pt_actual = pt_data(pt_i, pt_f, point_num, iter_num)
    pt_dev = abs(pt-pt_actual)
    plt.plot(pt_actual, pt_err, 'b.' , label="$p_t$ Error")
    plt.plot(pt_actual, pt_dev, 'r.' , label="$p_t$ Deviation")
    plt.xlabel('$p_t$')
    plt.ylabel('$\Delta p_t$') 
    plt.legend(loc='best',labelspacing=0.1 )
    plt.grid()
    plt.title('$p_t$ error vs $p_t$')
    plt.show()