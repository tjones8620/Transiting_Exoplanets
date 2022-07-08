#!/usr/bin/env python
# coding: utf-8
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



def cross_func(spectrum, spec_template):
    """
    Performs cross correlation on an inputed spectrum
    with respect to the spectral template and then 
    fits the data with a gaussian and outputs the mean
    of the gaussian
    """

    def correl_func(s, spec_template):
        """
        Defining a function that performs the scipy.signal.correlate() function
        """
        return correlate(s - s.mean(), spec_template - spec_template.mean(), mode = 'same')
    
    corr = correl_func(spectrum, spec_template)
    
    fit, var = curve_fit(gaussian, v, corr, p0 = [4000, 1000, 0, 1000])
    mean, sigma, b, a = fit
    return mean



def transit_flux(t, T_0, a, p, i, f_oot):
    
    
    #defining empty arrays to store x, y, z
    x, y, z = np.ones(len(t_tran)), np.ones(len(t_tran)), np.ones(len(t_tran))
    
    flux = np.ones(len(t_tran))     #defining array of length t to store flux values initialised at 1
    
    
    for j in range(len(t_tran)):
   
        phi = 2*np.pi/P * (t_tran[j] - T_0)  #converting each time to phase
                
        x[j], y[j] = a * np.sin(phi), a * np.cos(phi)  #calculating orbital coordinates and appending x,y lists
        
        
        z_t = np.sqrt(x[j]**2 + (np.cos(i) * y[j])**2)    #calculating z(t)
        z[j] = z_t
    
        k_0 = np.arccos((p**2 + z[j]**2 - 1)/(2 * p * z[j]))  #computing k0 and k1 values for each z
        k_1 = np.arccos((1 - p**2 + z[j]**2)/(2 * z[j]))

        
        #making conditions for flux
        if z[j] > 1 + p:
            flux[j] = 1
            
        if z[j] <= 1 - p:
            flux[j] = (1 - p**2)
            
        if 1 - p < z[j] <= 1 + p:  
            flux[j] = 1 - 1/np.pi*(p**2*k_0 + k_1 - np.sqrt((4*z[j]**2-(1+z[j]**2-p**2)**2)/4))

            
    return f_oot * flux  #returning flux array multiplied by out-of-transit flux




#defining sine function that can be used to fit sinusoidal data
def sin_fit(x, K, d):
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
    c = (((2*np.pi)/P)*(T_0+diff_in_time))

    return K * np.sin(b * x - c)+d


def residuals(y,f, x1, x2, x3, s=1):
    return (y - f(x1,x2,x3))/s


def gauss_error(A, x, dx, y, dy):
    return (A * np.sqrt((dx/x)**2 + (dy/y)**2))









