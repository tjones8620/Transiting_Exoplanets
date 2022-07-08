#!/usr/bin/env python
# coding: utf-8

from astro_modelling_functions import*

"""import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.optimize as opt
from scipy.signal import correlate
"""

 

#importing data for rv and transient
transit = np.load('JSAstroLab2022_transit_data.npz')
rv = np.load('JSAstroLab2022_rv_data.npz')


#defining mass radius and temperature of the star
M = 0.102
R = 0.137
Temp = 2800


rv_h = rv.files       #defining handles for rv

spectra = rv_h[3::2]  #defining handles for spectra and times
times = rv_h[4::2]

wl = rv[rv_h[0]]      #extracting wavelength, velocity and spectral template from the rv file
v = rv[rv_h[1]]
spec_template = rv[rv_h[2]]


from astro_modelling_functions import*

#test plot of one of the spectra against velocity superimposed on the spectral template
fig1, ax1 = plt.subplots(figsize = (12,8))
ax1.plot(wl, spec_template, 'k', label = "Spectral Template")
ax1.plot(wl, rv[spectra[10]], 'r.', color = 'cyan', label = "Recorded Spectrum")
#ax1.set_title("Spectrum")
ax1.set_xlabel("Wavelength (Angstrom)", fontsize = 20)
ax1.set_ylabel("Relative Flux")
ax1.legend()



#test plot for velocity shift of first spectrum with respect to the spectral template
spectrum = rv[spectra[0]]
corr = correlate(spectrum - spectrum.mean(), spec_template - spec_template.mean(), mode = 'same')

fit, var = curve_fit(gaussian, v, corr, p0 = [4000, 1000, 0, 1000])
mean, sigma, b, a = fit

fig2, ax2 = plt.subplots(figsize = (12,8))
ax2.plot(v, corr, 'k.', label = "Data", markersize = 8)
ax2.plot(v, gaussian(v, mean, sigma, b, a), 'cyan', label = "Gaussian Fit to Data", linewidth = 2)
ax2.set_xlim(-0.2e6, 0.2e6)
ax2.set_xlabel("Velocity Shift", fontsize = 15)
ax2.legend()
#print(v[np.argmax(corr)])





#empty list in which velocity shifts will be stored
shift = np.zeros(len(spectra))

#empty list in which times will be stored
time_list = np.zeros(len(times))

for j in range(len(spectra)):
    #loop over all spectra to find velocity shifts of every spectrum w.r.t the spectral template
    v_shift = cross_func(rv[spectra[j]], spec_template)
    shift[j] = v_shift
    
    #loop through rv data to store times using time handles
    time_list[j] = float(rv[times[j]]) #converting each size 1 numpy array into a float    




rv_starting_time = time_list.min()  #timestamp of first recorded spectrum in the rv data

time_list -= time_list.min() #bringing times back to zero for simplicity




#time array for which to plot the sinusoidal fit over
t = np.linspace(0.0,2.15,100)
#ax3.plot(t, sin_fit(t, K, c, d), color = 'red')

#plotting the individual v shifts against their respective timestamps
fig3, ax3 = plt.subplots(figsize = (12,8))
ax3.plot(time_list, shift, 'k.', markersize = 8)
ax3.set_xlabel("Time (days)", fontsize = 20)
ax3.set_ylabel(r"Velocity Shift $(ms^{-1})$", fontsize = 20)


# # Task 2 - Transit Model
# 
# The transit data was downloaded and the times and measured flux were stored in separate arrays. The period, $\textit{P}$, was stored in a separate variable. Firstly the measured flux of the star as a function of time was plotted in order to see the light curve of the star. It is clearly seen in this plot that the light curve drops dramatically, which corresponds to when the exoplanet passes in front of the star. 
# 
# Next a function is defined which plots a model of the transit. The scipy.optimize.curve_fit is used to produce a fit to the transit data using this transit function, outputting the fit parameters.
# 
# From the fit it is possible to determine the radius of the planet from the parameter $\rho$, using the expression:
# 
# $$ \rho = \frac{R_p}{R_{\star}} $$
# 
# using the fact that the radius of the star is given as $R_{\star} = 0.137R_{\odot}$. This is then multiplied by the ratio $R_{\odot}/R_{earth}$ to give the radius of the explanet in terms of solar masses. 
# 



# stores timestamps and respective flux in arrays
# stores period in a separate variable
t_tran, fl, P = transit['time'], transit['flux'], transit['P'] 

#defining the difference between the starting times of the transit and rv data
diff_in_time = t_tran.min() - rv_starting_time   

print(diff_in_time)

t_tran -= t_tran.min()  #bringing times back to start at zero for simplicity 




fig4, ax4 = plt.subplots(figsize = (12,8))
ax4.plot(t_tran, fl, 'k.')
ax4.set_xlabel("Time (days)", fontsize = 20)
ax4.set_ylabel("Relative Flux", fontsize = 20)
#ax4.set_title("Light curve")




fit3, var3 = curve_fit(transit_flux, t_tran, fl, p0 = [0.14, 15, 0.1, np.pi/2+0.001, 1] )  #fitting the transit function
T_0, a, p, i, f_oot = fit3  #transit parameters
print(fit3)

fit3_errors = np.sqrt(np.diag(var3))  #errors in the transit parameters

fig5, ax5 = plt.subplots(figsize = (12,8))
ax5.plot(t_tran, fl, 'k.', label = "Transit Data")
ax5.plot(t_tran, transit_flux(t_tran, T_0, a, p, i, f_oot), 'cyan', linewidth = 2, label = "Fit to Transit Data")
#ax5.set_title("Transit Model")
ax5.set_xlabel("Time (days)", fontsize = 20)
ax5.set_ylabel("Relative Flux", fontsize = 20)
ax5.legend()




#extracting the error in the parameters for the  from the fit
dT_0 = fit3_errors[0]
da = fit3_errors[1]
dp = fit3_errors[2]
di = fit3_errors[3]
dsini = abs(np.cos(i) * di)   #error in sin(i) using gauss' error law
df_oot =fit3_errors[4]




r = R * p   #finding the radius of the planet in solar radii
dr = R * dp  #error in radius of the planet in solar radii
Earth_Sun_radius = 109.076  #ratio of sun radius to earth radius

r_p = r * Earth_Sun_radius  #writing the radius of the planet in earth radii
dr_p = dr * Earth_Sun_radius #error in the radius of the planet in earth radii
print(f"Radius of expoplanet = {r_p:.2f} +- {dr_p:.2f} Earth Radii")

r_earth_m = 6.371e+6 #earth radius in meters
r_m = r_p * r_earth_m  #writing the radius of the planet in metres
dr_m = dr_p * r_earth_m #error in radius of the planet in metres


# # Task 3 
#         
# 



G = 6.67408e-11  #gravitational constant




folded_times = np.zeros(len(time_list))

folded_times = np.where(time_list < P, time_list, time_list - P)    #folding times to cover just one period




fit4, var4 = curve_fit(sin_fit, folded_times, shift, p0 = [-30, 5000])   #fitting for the rv curve
K, d = fit4  #rv curve fit parameters

#
fig6, ax6 = plt.subplots(figsize = (12,8))
t = np.linspace(0, P, 100)

ax6.plot(t, sin_fit(t, K, d), color = 'cyan', linewidth = 3, label = "Fit of RV Data")
ax6.plot(folded_times, shift, 'k.', markersize = 8, label = "RV Data")
#ax6.plot(t, sin_fit(t, K, ((2*np.pi)/P)*(T_0+diff_in_time)-np.pi, d), color = 'red')
#ax6.set_title("Phase folded plot of Velocity Shift against Time")
ax6.set_xlabel("Time (days)", fontsize = 20)
ax6.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 20)
ax6.legend()
print(fit4)
print(np.sqrt(np.diag(var4)))




#print(residuals(shift, sin_fit(time_list, K,b,c,d)))
resid = residuals(shift, sin_fit, folded_times, K, d)

fig7, (ax, ax1) = plt.subplots(nrows = 2, figsize = (12,8))
ax.plot(t, sin_fit(t, K, d), color = 'cyan', linewidth = 3, label = "Fit of RV Data")
ax.plot(folded_times, shift, 'k.', markersize = 8, label = "RV Data")
ax.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 20)

ax1.plot(folded_times, resid, 'k.', markersize = 8)
ax1.axhline(0, 0, 1.2, color = 'cyan', linewidth =3)
ax1.set_xlabel("Time (days)", fontsize = 20)

RMS = np.std(resid)
print(RMS)




#fitting for rv data taking into account the error in the data
fit4, var4 = curve_fit(sin_fit, folded_times, shift, p0 = [-28, 5000], sigma = np.ones(len(shift))*RMS, absolute_sigma = True)
K, d = fit4  #new fitting parameters for rv curve 

print(K, d)
print(np.sqrt(np.diag(var4)))

dK = np.sqrt(np.diag(var4)[0])  #error in the amplitude of the rv curve fit

t = np.linspace(0, P, 100)

fig8, ax = plt.subplots(figsize = (12,8))
ax.plot(t, sin_fit(t, K, d), color = 'cyan', linewidth = 3, label = "Fit to RV Data")
ax.errorbar(folded_times, shift, fmt = 'k.', yerr = RMS*np.ones(len(shift)), capsize = 3, markersize = 8 , label = "RV Data")
#ax7.set_title("Phase folded plot of Velocity Shift against Time")
ax.set_xlabel("Time (days)", fontsize = 20)
ax.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 20)
ax.legend()



shift_nosystemic = shift - d  #taking away systemic velocity to see the true values of velocity shift




fig9, ax8 = plt.subplots(figsize = (12,8))
ax8.plot(t, sin_fit(t, K, 0), color = 'cyan', linewidth =3 , label = "Fit to RV Data")
ax8.errorbar(folded_times, shift_nosystemic, yerr = RMS, fmt = 'k.', capsize = 3, markersize = 8, label = "RV Data")
#ax8.set_title("Phase folded plot of Velocity Shift against Time")
ax8.set_xlabel("Time (days)", fontsize = 20)
ax8.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 20)
ax8.legend()




K_abs = abs(K)  #taking absolute value of amplitude 

Msun_kg = 1.99e30 #mass of sun in kg
M_kg = M *  Msun_kg #mass of this star in kg

m_kg = ((((M_kg**(2/3)) * K_abs)/np.sin(i)) * ((P*24*60*60)/(2 * np.pi *G))**(1/3))  #mass of the planet in kg
dm_kg = gauss_error(m_kg ,np.sin(i),  dsini, K_abs, dK) #error in planet mass in kg from gauss' error law

M_earth = 5.972e24  #mass of earth in kg
m = m_kg/M_earth #mass of planet in earth masses
dm = gauss_error(m ,np.sin(i),  dsini, K_abs, dK)  #error in planet mass in earth masses

min_m = m*np.sin(i)
dmin_m = gauss_error(min_m, m, dm, np.sin(i), dsini)

print(f"Mass of exoplanet = {m:.5f} +- {dm:.5f} Earth Masses")
print(f"Minimum mass of exoplanet = {min_m:.5f} +- {dmin_m:.5f} Earth Masses")

solar_au = 0.0046524726 #solar radius in au
R_star_au = solar_au * R  #radius of this star in au
semi_major = R_star_au * a  #semi major axis of this system in au
dsemi_major = R_star_au * da  #error in semi major axis in au
print(f"Semi major axis of planet-star system ={semi_major:.4f} +- {dsemi_major:.4f}au")


vol_p = (4/3 * np.pi * (r_m)**3)  #volume of the planet in m^3
dvol_p = (vol_p * 3 * 4/3 * np.pi * dr_m)/r_m #error in planet volume in m^3

density = m_kg/vol_p  #density of planet in kg/m^3
d_density = gauss_error(density, vol_p, dvol_p, m_kg, dm_kg) #error in density of planet in kg/m^3

earth_density = 5520 #in kg

dens_p = density/earth_density #density of the planet in earth densities
ddens_p = d_density/earth_density # error in density of planet in earth densities


print(f"Density of exoplanet = {density:.2f} +- {d_density:.2f} kg/m^3")
print(f"Density of exoplanet = {dens_p:.2f} +- {ddens_p:.2f} earth densities")









t1 = t_tran + diff_in_time  #starting time of transit data w.r.t rv data

fig10, (ax9, ax10) = plt.subplots(nrows = 2, figsize = (12, 8))

ax9.plot(t_tran, transit_flux(t_tran, T_0, a, p, i, f_oot), 'cyan')
ax9.axvline(x = T_0, color = 'k', linestyle = '--')
ax9.set_ylabel("Relative Flux", fontsize = 15)

ax10.plot(t1, sin_fit(t1, K, 0), color = 'cyan')
ax10.axvline(x = T_0 + diff_in_time, color = 'k', linestyle = '--')
ax10.axhline(0, 0, 1.1, color = 'k', linestyle = '--')
ax10.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 15)
ax10.set_xlabel("Time (days)", fontsize = 15)

fig9.subplots_adjust(hspace = 0)









fig1.savefig("test_spectrum.jpg")
fig2.savefig("v_correlate.jpg")
fig3.savefig("rv_data.jpg")
fig4.savefig("transit_data.jpg")
fig5.savefig("transit_fit.jpg")
fig6.savefig("rv_fit.jpg")
fig8.savefig("rvfit_errors.jpg")
fig9.savefig("rvfit_nosystemic.jpg")
fig10.savefig("rv_vs_transit.jpg")


