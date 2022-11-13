from exoplanets_functions import *
import astropy.units as u
import astropy.constants as const
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

########### Importing rv and transit data ###########
transit = np.load('JSAstroLab2022_transit_data.npz')
rv = np.load('JSAstroLab2022_rv_data.npz')

############ Defining mass, radius and temperature of the star ########
M = 0.102
R = 0.137
Temp = 2800
G = const.G.value #gravitational constant

######################################## PART 1: RV DATA ############################################
########################### Finding and plotting radial velocity data ###############################
#####################################################################################################

rv_h = rv.files       #defining handles for rv
spectra = rv_h[3::2]  #defining handles for spectra and times
times = rv_h[4::2]

# Extracting wavelength, velocity and spectral template from the rv file
wl = rv[rv_h[0]]      
v = rv[rv_h[1]]
spec_template = rv[rv_h[2]]

# Test plot of one of the spectra against velocity superimposed on the spectral template
fig1, ax1 = plt.subplots(figsize = (12,8))
ax1.plot(wl, spec_template, 'k', label = "Spectral Template")
ax1.plot(wl, rv[spectra[0]], 'c.', label = "Recorded Spectrum", ms=2)
ax1.set_xlabel("Wavelength (Angstrom)", fontsize = 20)
ax1.set_ylabel("Relative Flux", fontsize = 20)
ax1.legend()

#test plot for velocity shift of first spectrum with respect to the spectral template
spectrum = rv[spectra[0]]  #selecting the first spectrum in the dataset to plot

#correlation method with the spectrum and the spectral template
corr = correlate(spectrum - spectrum.mean(), spec_template - spec_template.mean(), mode = 'same') 

#fitting a gaussian function to the velocity shift using curve_fit
fit, var = curve_fit(gaussian, v, corr, p0 = [4000, 1000, 0, 1000])
mean, sigma, b, a = fit #defining the fitting parameters

print("Section 1: Radial Velocity Data")
print(f"Mean value of velocity for this measured spectrum is {mean:.3f} m/s and error in the mean is {np.sqrt(np.diag(var)[0]):.3f} m/s")

#plotting this cross correlation of the spectrum to the template with the fitted gaussian function
fig2, ax2 = plt.subplots(figsize = (12,8))
ax2.plot(v, corr, 'k.', label = "Data", markersize = 8)
ax2.plot(v, gaussian(v, mean, sigma, b, a), 'cyan', label = "Gaussian Fit to Data", linewidth = 2)
ax2.set_xlim(-0.2e6, 0.2e6)
ax2.set_xlabel("Velocity Shift", fontsize = 15)
ax2.legend()


shift = np.zeros(len(spectra)) #empty list in which velocity shifts will be stored
time_list = np.zeros(len(times)) #empty list in which times will be stored

#loop over all spectra to find shift of every spectrum w.r.t the spectral template and store times in an array
for j in range(len(spectra)):   
    v_shift = cross_func(rv[spectra[j]], spec_template, v= v)
    shift[j] = v_shift 
    time_list[j] = float(rv[times[j]]) #converting each size 1 numpy array into a float   
    
    
rv_starting_time = time_list.min()  #timestamp of first recorded spectrum in the rv data
time_list -= time_list.min() #bringing times back to zero for simplicity

#time array for which to plot the sinusoidal fit over
t = np.linspace(0.0,2.15,100)

#plotting the individual radial velocities against their respective timestamps
fig3, ax3 = plt.subplots(figsize = (12,8))
ax3.plot(time_list, shift, 'k.', markersize = 8)
ax3.set_xlabel("Time (days)", fontsize = 20)
ax3.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 20)

input("Press Enter to continue to Part 2: Transit Data")

####################################### PART 2: TRANSIT DATA ########################################
################################# Finding and plotting transit data #################################
#####################################################################################################

# stores timestamps and respective flux in arrays and period in a variable P
t_tran, fl, P = transit['time'], transit['flux'], transit['P'] 

# Defining the difference between the starting times of the transit and rv data
diff_in_time = t_tran.min() - rv_starting_time   

print(f"Section 2: Transit Data")
print(f"Difference between the starting time of the transit and rv data is {diff_in_time:.3f} days")

t_tran -= t_tran.min()  #bringing times back to start at zero for simplicity 

#plotting the recorded flux as a function of time
fig4, ax4 = plt.subplots(figsize = (12,8))
ax4.plot(t_tran, fl, 'k.')
ax4.set_xlabel("Time (days)", fontsize = 20)
ax4.set_ylabel("Relative Flux", fontsize = 20)

print(P)
#fitting the transit data using the transit model
bounds = ([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, P-0.000001], [np.inf, np.inf, np.inf, np.inf, np.inf, P+0.000001])
fit3, var3 = curve_fit(transit_flux, t_tran, fl, p0 = [0.14, 15, 0.1, np.pi/2+0.001, 1, P] , bounds=bounds)  #fitting the transit function
T_0, a, p, i, f_oot, _ = fit3  #transit parameters
fit3_errors = np.sqrt(np.diag(var3)[:-1])  #errors in the transit parameters

print(f"Fitted transit parameters T_0, a, p, i, f_oot are {fit3[:-1]} with errors {fit3_errors[:-1]}")


#plotting the fit to the transit model with the recorded transit data
fig5, ax5 = plt.subplots(figsize = (12,8))
ax5.plot(t_tran, fl, 'k.', label = "Transit Data")
ax5.plot(t_tran, transit_flux(t_tran, T_0, a, p, i, f_oot, P), 'cyan', linewidth = 2, label = "Fit to Transit Data")
ax5.set_xlabel("Time (days)", fontsize = 20)
ax5.set_ylabel("Relative Flux", fontsize = 20)
ax5.legend()


# Extracting the errors in the fit parameters
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

r_earth_m = const.R_earth.to('m')  #radius of earth in metres
r_m = r_p * r_earth_m  #writing the radius of the planet in metres
dr_m = dr_p * r_earth_m #error in radius of the planet in metres

input("Press Enter to continue to Part 3: Transit Data")
############################################################################################################
#################################### PART 3: TRANSIT AND RV DATA ###########################################
############################################################################################################

folded_times = np.zeros(len(time_list))

folded_times = np.where(time_list < P, time_list, time_list - P)    #folding times to cover just one period

#fitting a sine function to the rv curve using
bounds = ([-np.inf, -np.inf, T_0+diff_in_time - 0.000001, P-0.000001], [np.inf, np.inf, T_0+diff_in_time + 0.000001, P+0.000001])
fit4, var4 = curve_fit(sin_fit, folded_times, shift, p0 = [-30, 5000, T_0+diff_in_time, P], bounds=bounds)  #fitting the sine function
K, d, _, _ = fit4  #rv curve fit parameters

print(f"Amplitude of fit is {fit4[0]:.3f} +- {np.sqrt(np.diag(var4))[0]:.3f} ")

t = np.linspace(0, P, 100) #array for which to plot the rv fit against

#plotting the fit of the rv curve with the rv data
fig6, ax6 = plt.subplots(figsize = (12,8))
ax6.plot(t, sin_fit(t, K, d, T_0+diff_in_time, P), color = 'cyan', linewidth = 3, label = "Fit of RV Data")
ax6.plot(folded_times, shift, 'k.', markersize = 8, label = "RV Data")
ax6.set_xlabel("Time (days)", fontsize = 20)
ax6.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 20)
ax6.legend()


resid = residuals(shift, sin_fit, folded_times, K, d, T_0+diff_in_time, P) #finding the residuals of the data w.r.t the best fit


#plotting the rv curve and residuals as functions of time
fig7, (ax, ax1) = plt.subplots(nrows = 2, figsize = (12,8))
ax.plot(t, sin_fit(t, K, d, T_0+diff_in_time, P), color = 'cyan', linewidth = 3, label = "Fit of RV Data")
ax.plot(folded_times, shift, 'k.', markersize = 8, label = "RV Data")
ax.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 20)
ax1.plot(folded_times, resid, 'k.', markersize = 8)
ax1.axhline(0, 0, 1.2, color = 'cyan', linewidth =3)
ax1.set_xlabel("Time (days)", fontsize = 20)

RMS = np.std(resid) #taking the standard deviation of the residuals to be the common error of the data
print(RMS)

#fitting for rv data taking into account the error in the data
bounds = ([-np.inf, -np.inf, T_0+diff_in_time - 0.000001, P-0.000001], [np.inf, np.inf, T_0+diff_in_time + 0.000001, P+0.000001])
fit5, var5 = curve_fit(sin_fit, folded_times, shift, p0 = [-28, 5000, T_0+diff_in_time, P], sigma = np.ones(len(shift))*RMS, absolute_sigma = True)
K, d, _, _ = fit5  #new fitting parameters for rv curve 

print(f"RV curve fitted parameters are K = {K:.3f}+-{np.sqrt(np.diag(var5)[:-2])[0]:.3f}, d = {d:.3f}+-{np.sqrt(np.diag(var5)[:-2])[1]:.3f}")

dK = np.sqrt(np.diag(var5)[0])  #error in the amplitude of the rv curve fit


#plotting rv curve and best fit including the errors in the data
fig8, ax = plt.subplots(figsize = (12,8))
ax.plot(t, sin_fit(t, K, d, T_0+diff_in_time, P), color = 'cyan', linewidth = 3, label = "Fit to RV Data")
ax.errorbar(folded_times, shift, fmt = 'k.', yerr = RMS*np.ones(len(shift)), capsize = 3, markersize = 8 , label = "RV Data")
ax.set_xlabel("Time (days)", fontsize = 20)
ax.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 20)
ax.legend()


shift_nosystemic = shift - d  #taking away systemic velocity to see the true values of velocity shift

#plotting the fit of the rv curve with the systemic velocity subtracted
fig9, ax8 = plt.subplots(figsize = (12,8))
ax8.plot(t, sin_fit(t, K, 0, T_0+diff_in_time, P), color = 'cyan', linewidth =3 , label = "Fit to RV Data")
ax8.errorbar(folded_times, shift_nosystemic, yerr = RMS, fmt = 'k.', capsize = 3, markersize = 8, label = "RV Data")
ax8.set_xlabel("Time (days)", fontsize = 20)
ax8.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 20)
ax8.legend()


#calculating the mass of the planet
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


input("Press enter to continue to Part 4: System Characteristics")
############################################################################################################
###################################  PART 4: SYSTEM CHARACTERISTICS ########################################
############################################################################################################

#finding the semimajor axis of the system
au = 149597871000 #1au in metres
solar_au = 0.0046524726 #solar radius in au
R_star_au = solar_au * R  #radius of this star in au
semi_major = R_star_au * a  #semi major axis of this system in au
dsemi_major = R_star_au * da  #error in semi major axis in au
print(f"Semi major axis of planet-star system ={semi_major:.4f} +- {dsemi_major:.4f}au")

#calculating the density of the planet

vol_p = (4/3 * np.pi * (r_m)**3)  #volume of the planet in m^3
dvol_p = (vol_p * 3 * 4/3 * np.pi * dr_m)/r_m #error in planet volume in m^3

density = m_kg/vol_p  #density of planet in kg/m^3
d_density = gauss_error(density, vol_p, dvol_p, m_kg, dm_kg) #error in density of planet in kg/m^3

earth_density = 5520 #in kg

dens_p = density/earth_density #density of the planet in earth densities
ddens_p = d_density/earth_density # error in density of planet in earth densities


print(f"Density of exoplanet = {density:.2f} +- {d_density:.2f} kg/m^3")
print(f"Density of exoplanet = {dens_p:.2f} +- {ddens_p:.2f} earth densities")


#Calculating the flux at the surface of the exoplanet
R_star_m = R_star_au * au #radius of the star in m
stefboltz = 5.6704e-8 #stefan boltzmann constant
surface_flux = stefboltz * Temp**4 #surface flux of the star
luminosity = surface_flux * 4*np.pi*(R_star_m)**2 #luminosity of the star in W
semi_major_m = semi_major*au #semi major axis in m
surface_flux = luminosity/(4*np.pi*semi_major_m**2) #flux at surface of exoplanet in W/m^2
flux_at_earth = 1361 #flux at earths surface


print(surface_flux, "W/m^2")
print(surface_flux/flux_at_earth)


t1 = t_tran + diff_in_time  #starting time of transit data w.r.t rv data


#plotting transit curve and the part of the rv curve that corresponds to the same time period
fig10, (ax9, ax10) = plt.subplots(nrows = 2, figsize = (12, 8))
ax9.plot(t_tran, transit_flux(t_tran, T_0, a, p, i, f_oot, P), 'cyan')
ax9.axvline(x = T_0, color = 'k', linestyle = '--')
ax9.set_ylabel("Relative Flux", fontsize = 15)
ax10.plot(t1, sin_fit(t1, K, 0, T_0+diff_in_time, P), color = 'cyan')
ax10.axvline(x = T_0 + diff_in_time, color = 'k', linestyle = '--')
ax10.axhline(0, 0, 1.1, color = 'k', linestyle = '--')
ax10.set_ylabel(r"Radial Velocity $(ms^{-1})$", fontsize = 15)
ax10.set_xlabel("Time (days)", fontsize = 15)
fig9.subplots_adjust(hspace = 0)