from mu import *
from scipy import optimize

xhit, yhit, zhit = np.asarray(muon_path(plot = False))
rhit = np.sqrt(xhit*xhit + yhit*yhit)
phihit = np.arctan2(yhit, xhit)     

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

err_layers = [err_1, err_2, err_3, err_4]
det_layers = [[], r_1, r_2, r_3, r_4]

def smear_func(initial, final, err):
    i = len(initial)
    f = i + len(final)
    r_smear = rhit[i:f]
    phi_smear = np.random.normal(phihit[i:f], err[0]/rhit[i:f])
    x_smear = r_smear * np.cos(phi_smear)
    y_smear = r_smear * np.sin(phi_smear)
    z_smear = np.random.normal(zhit[i:f], err[1])
    return x_smear, y_smear, z_smear

    
def smear_data():
    x_smear, y_smear, z_smear = [],[],[]
    err_layers = [err_1, err_2, err_3, err_4]
    det_layers = [[], r_1, r_2, r_3, r_4]
    for i in range(len(det_layers)-1):
        x_smear = np.append(x_smear, smear_func(det_layers[i], det_layers[i+1], err_layers[i])[0])
        y_smear = np.append(y_smear, smear_func(det_layers[i], det_layers[i+1], err_layers[i])[1])
        z_smear = np.append(z_smear, smear_func(det_layers[i], det_layers[i+1], err_layers[i])[2])
    return x_smear, y_smear, z_smear

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
    residu   = np.sum((Ri - R)**2)
    return xc, yc, R, residu
 
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

x_data, y_data, z_data = smear_data()
xc, yc, R, resid = leastsq_circle(x_data, y_data)
plot_data_circle(x_data, y_data, xc, yc, R)
print xc, yc, R