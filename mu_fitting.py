from mu import *
from scipy import optimize

r_pixel = [0.04,0.07,0.11]
r_tib=[0.26,0.32,0.43,0.62]
r_tob=[0.71,0.79,0.88,0.97,1.07,1.14]
r_tracker=np.append(np.append(r_pixel,r_tib), r_tob)


r_mb1= [4.20,4.22,4.24,4.26,4.42,4.44,4.46,4.48]
r_mb2= [5.00,5.02,5.04,5.06,5.24,5.26,5.28,5.30]
r_mb3= [6.08,6.10,6.12,6.14,6.32,6.34,6.36,6.38]
r_mb4= [7.10,7.12,7.14,7.16,7.34,7.36,7.38,7.40]

r_drift=np.append(np.append(r_mb1,r_mb2),np.append(r_mb3,r_mb4))

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
err_4_t=100.*10**(-6) # MB transverse error, units: meters
err_4_z=150.*10**(-6) # MB z error, units: meters
err_4 = [err_4_t, err_4_z]


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
    f = i + final
    r_smear = rhit[i:f]
    phi_smear = np.random.normal(phihit[i:f], err[0]/rhit[i:f])
    x_smear = r_smear * np.cos(phi_smear)
    y_smear = r_smear * np.sin(phi_smear)
    z_smear = np.random.normal(zhit[i:f], err[1])
    return f, x_smear, y_smear, z_smear
    

#Calculates the actual hits of the particle in case its initial position is
#beyond the innermost tracker.
def determine_start(rhit):
    r_pixel_new = []
    r_tib_new = []
    r_tob_new = []
    r_drift_new = []
    for j in range(len(r_pixel)):
        if rhit[0]<(r_pixel[j]+0.01):
            r_pixel_new = np.append(r_pixel_new, r_pixel[j])
    if len(r_pixel_new)>0:
        return r_pixel_new, r_tib, r_tob, r_drift
    else:
        for j in range(len(r_tib)):
            if rhit[0]<(r_tib[j]+0.01):
                r_tib_new = np.append(r_tib_new, r_tib[j])
        if len(r_tib_new)>0:
            return r_pixel_new, r_tib_new, r_tob, r_drift   
        else:
            for j in range(len(r_tob)):
                if rhit[0]<(r_tob[j]+0.01):
                    r_tob_new = np.append(r_tob_new, r_tob[j])
            if len(r_tob_new)>0:
                return r_pixel_new, r_tib_new, r_tob_new, r_drift
            else:
                for j in range(len(r_drift)):
                    if rhit[0]<(r_drift[j]+0.01):
                        r_drift_new = np.append(r_drift_new, r_drift[j])
                return r_pixel_new, r_tib_new, r_tob_new, r_drift_new


#Smears each data point in the appropriate way (data from multiple trackers).    
def smear_data(rhit, phihit, zhit, tube):
    x_smear, y_smear, z_smear = [],[],[]
    err_layers = [err_1, err_2, err_3, err_4]
    det_layers = [[], r_pixel, r_tib, r_tob, r_drift]
    if rhit[0]>r_pixel[0]: 
        temporary = determine_start(rhit)
        det_layers = [[], temporary[0], temporary[1], temporary[2], temporary[3]]
        r_tracker2 = np.concatenate((det_layers[0], det_layers[1], det_layers[2], det_layers[3]))
    temp = [0]
    if tube == 1:
        for j in range(len(det_layers)-2):
            temp = smear_func(temp[0], len(det_layers[j+1]), err_layers[j], rhit, phihit, zhit)
            x_smear = np.append(x_smear, temp[1])
            y_smear = np.append(y_smear, temp[2])
            z_smear = np.append(z_smear, temp[3])
    if tube == 2:
        if rhit[0] > r_pixel[0]: temp = smear_func(len(r_tracker2), len(det_layers[4]), err_layers[3], rhit, phihit, zhit)
        else: temp = smear_func(len(r_tracker), len(det_layers[4]), err_layers[3], rhit, phihit, zhit)
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
def pt_calc(rhit, phihit, zhit,tube):
    x_data, y_data, z_data = smear_data(rhit, phihit, zhit, tube)
    if len(x_data)>2: xc, yc, Rc = leastsq_circle(x_data, y_data)
    if len(x_data)<=2: Rc=0 #Cannot fit circle with 0,1, or 2 points.
    if tube == 1:
        pt = 0.3*B*Rc
    if tube == 2:
        pt = 0.3*B2*Rc
    return pt

#Iterates calculation of transverse momentum over a given amount of times,
#and calculates a mean and a standard error for the estimated transverse
#momentum for this particle.
def pt_datapoint(pt, iter_num):
    rhit, phihit, zhit = gen_data(pt=pt)
    pt_data,pt_data2 = [],[]
    for i in range(iter_num):
        pt_data = np.append(pt_data, pt_calc(rhit, phihit, zhit,1))
        pt_data2 = np.append(pt_data2, pt_calc(rhit, phihit, zhit,2))
    pt_mean = np.mean(pt_data)
    pt_err = (np.std(pt_data))/(np.sqrt(iter_num))
    pt_mean2 = np.mean(pt_data2)
    pt_err2 = (np.std(pt_data2))/(np.sqrt(iter_num))
    return pt_mean, pt_err, pt_mean2, pt_err2


#Iterates previous function over a given range of transverse momenta.
def pt_data(pt_i, pt_f, point_num, iter_num):
    pt, pt_actual, pt_err, pt2, pt_err2 = [], [], [], [], []
    stepsize = float((pt_f - pt_i))/point_num
    for i in range(point_num):
        temp1, temp2, temp3, temp4 = pt_datapoint(iter_num = iter_num, pt = pt_i + i*stepsize)
        pt = np.append(pt, temp1)
        pt_err = np.append(pt_err, temp2)
        pt2 = np.append(pt2, temp3)
        pt_err2 = np.append(pt_err2, temp4)
        pt_actual = np.append(pt_actual, pt_i +i*stepsize)
    return pt, pt_err, pt2, pt_err2, pt_actual
    
#Plots the standard error in the transverse momentum and its deviation from
#the actual value against the transverse momentum.
def plot_data(pt_i=3, pt_f=20, point_num=100, iter_num=50): #Does not work if muon trapped in tracker.
    pt, pt_err, pt2, pt_err2, pt_actual = pt_data(pt_i, pt_f, point_num, iter_num)
    if all(pt) == 0: pt_dev = [0]*len(pt)
    else: pt_dev = abs(pt-pt_actual)
    pt2_dev = abs(pt2-pt_actual)
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(pt_actual, pt_err, 'b.' , label="$p_t$ Error")
    axarr[0].plot(pt_actual, pt_dev, 'r.' , label="$p_t$ Deviation")
    axarr[1].plot(pt_actual, pt_err2, 'b.' , label="$p_t$ Error")
    axarr[1].plot(pt_actual, pt2_dev, 'r.' , label="$p_t$ Deviation")
    plt.xlabel('$p_t$')
    plt.ylabel('$\Delta p_t$')
    plt.legend(loc='best',labelspacing=0.1,numpoints=1)
    plt.grid()
    plt.title('$p_t$ error vs $p_t$')
    plt.show()
    