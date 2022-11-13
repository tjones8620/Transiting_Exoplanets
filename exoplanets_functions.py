import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.optimize as opt
from scipy.signal import correlate

def gaussian(x, mean, sigma, b, a):
    """
    Gaussian function 
    """
    return a / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * (x - mean)**2 / sigma**2) + b

def cross_func(spectrum, spec_template, **kwargs):
    """
    Performs cross correlation on an inputed spectrum
    with respect to the spectral template and then 
    fits the data with a gaussian and outputs the mean
    of the gaussian
    """

    v = kwargs['v']

    def correl_func(s, spec_template):
        """
        Defining a function that performs the scipy.signal.correlate() function
        """
        return correlate(s - s.mean(), spec_template - spec_template.mean(), mode = 'same')
    
    corr = correl_func(spectrum, spec_template)
    
    fit, var = curve_fit(gaussian, v, corr, p0 = [4000, 1000, 0, 1000])
    mean, sigma, b, a = fit
    return mean



def transit_flux(t, T_0, a, p, i, f_oot, P):
    
    """
    - Takes in times t and parameters T_0, a, p, i, f_oot
    - Times t are changed into phase phi
    - x,y coordinates are defined
    - z(t) is defined in terms of x, y and orbital inclination i
    - k_0 and k_1 values are computed for each k
    - Conditonal statements determine the value that is added to the flux array
    - Function outputs this flux times f_oot
    """    



    flux = np.ones(len(t))     #defining array of length t to store flux values initialised at 1
    
    for j in range(len(t)):
   
        phi = 2*np.pi/P * (t[j] - T_0)  #converting each time to phase
                
        x, y = a * np.sin(phi), a * np.cos(phi)  #calculating orbital coordinates and appending x,y lists
         
        z_t = np.sqrt(x**2 + (np.cos(i) * y)**2)    #calculating z(t)
        z = z_t
    
        k_0 = np.arccos((p**2 + z**2 - 1)/(2 * p * z))  #computing k0 and k1 values for each z
        k_1 = np.arccos((1 - p**2 + z**2)/(2 * z))

        #making conditions for flux
        if z > 1 + p:
            flux[j] = 1
            
        if z <= 1 - p:
            flux[j] = (1 - p**2)
            
        if 1 - p < z <= 1 + p:  
            flux[j] = 1 - 1/np.pi*(p**2*k_0 + k_1 - np.sqrt((4*z**2-(1+z**2-p**2)**2)/4))

            
    return f_oot * flux  #returning flux array multiplied by out-of-transit flux

#defining sine function that can be used to fit sinusoidal data
def sin_fit(x, K, d, t, P):
    """
    Function for plotting sinusoidal curves
    
    b is given by 2pi/period; period is given
    at the start of this question in this question 
    so use this rather than fitting b as a parameter
    
    c is the phase difference. Can make use of 
    central transit time parameter, T_0 found with 
    transit curve fit. This has to be alterred 
    slightly to take into account the difference 
    in the t = 0 time between the transit data and 
    the radial velocity data.
    """
    
    b = (2 * np.pi / P) 
    c = (((2*np.pi)/P)*(t))

    return K * np.sin(b * x - c)+d


def residuals(y,f, x1, x2, x3, x4, x5, s=1):
    """
    producing the residuals of data
    with respect to the best fit
    """
    return (y - f(x1,x2,x3,x4,x5))/s

def gauss_error(A, x, dx, y, dy):
    """
    Gauss' cumulitive error for 
    multiplying errors
    """
    return (A * np.sqrt((dx/x)**2 + (dy/y)**2))