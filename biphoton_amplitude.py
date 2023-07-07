#-------------------------------------------------------------------------------
# Name:        SPDC - Phase matching
# Purpose:     PEEC
#
# Author:      TECM
#
# Created:     01/07/2023
# Copyright:   (c) TECM 2023
#-------------------------------------------------------------------------------

'''
REFERENCES:

    [1] Vittorio Giovannetti et al. "Extended phase-matching conditions for
    improved entanglement generation",


    [2] C. Chen. "Generation and characterization of spectrally
    factorable biphotons", MSc Thesis.

    [3] Vittorio Giovannetti et al. "Generating Entangled Two-Photon
    States with Coincident Frequencies"

'''




#===============================================================================
# PACKAGES
#===============================================================================

from pylab import *

#===============================================================================
#///////////////////////////////////////////////////////////////////////////////
#===============================================================================


#===============================================================================
# Functions
#===============================================================================

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
def refractive_index_BBO(omega):
    '''
    Refractive index of BBO

    Input (SI units):
        omega: angular frequency (float)
    '''

    wave0 = 2*pi*3e8/omega # convert to wavelength

    n = sqrt(2.7405 + 0.0184/(wave0**2-0.0179) - 0.0155*wave0**2)

    return n
#///////////////////////////////


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
def Gaussian_spec(w_array, wp, Omega_p):
    '''
    Guassian spectrum

    Inputs (SI units):
        w_array: angular frequency array (float or array)
        wp:      pump central frequency  (float)
        Omega_p: pump spectral width     (float)
    '''
    return exp(-(w_array-wp)**2/(Omega_p**2))
#///////////////////////////////


#///////////////////////////////
def alpha_func(ws, wi, ns, ni, qsi_p):
    '''
    Amplitude biphoton

    Inputs (SI units):
        ws: angular frequency of signal (float)
        wi: angular frequency of idler (float)
        ns: refractive index of signal (float)
        ni: refractive index of idler (float)
        qsi_p: pump spectrum (float)
    '''

    res = qsi_p #* sqrt(ws*wi)/(ns*ni)
    return res
#///////////////////////////////


#///////////////////////////////
def Phi_L_func(L, wp, ws, wi, gamma_s, gamma_i, theta):
    '''
    Phase Matching Function

    Inputs (SI units):
        L: Crystal length (float)
        wp: angular frequency of pump (float)
        ws: angular frequency of signal (float)
        gamma_s:
        gamma_i:
        theta:   angle (float)

    '''

    dk = (ws-wp/2)*gamma_s*cos(theta) + (wi-wp/2)*gamma_i*sin(theta)

    res = sinc(dk*L/(2*pi))  # note: this extra pi is because of the numpy's definition of sync

    return res
#///////////////////////////////


#///////////////////////////////
def wavenumber(omega, n):
    '''
    Wavenumber calculation

    Input:
        omega: angular frequency (float)
        n    : refractive index  (float)
        '''
    return omega * n/3e8
#///////////////////////////////


#===============================================================================
#///////////////////////////////////////////////////////////////////////////////
#===============================================================================



#===============================================================================
# ARRAYS
#===============================================================================



#----------------------------
# Parameters (SI units)
#---------------------------

# array size
N = 300

# wavelengths
wavelength_1 = 705.357e-9
wavelength_2 = 897.727e-9
wavelength_c = 900e-9

# angular frequencies
w1 = 2*pi*3e8/wavelength_1
w2 = 2*pi*3e8/wavelength_2
wc = 2*pi*3e8/wavelength_c

## arrays

# signal array
ws_array = linspace(w2, w1, N)

# idler array
wi_array = linspace(w2, w1, N)

# pump array
wp_array = linspace(w2, w1, N)

#===============================================================================
#///////////////////////////////////////////////////////////////////////////////
#===============================================================================



#===============================================================================
# SIMULATION 1
#===============================================================================


#-----------------
# Parameters (SI units)
#-----------------
L       = 1e-2
theta   = pi/20
Omega_p = 6e13
gamma   = 8e-5 * 1e-12/1e-6

# pump central angular frequency spectrum
wp      = 2*pi*3e8/395e-9    # lambda_0 = 395e-9
#-----------------


#------------------------------------------------
# SIMULATION
#------------------------------------------------

# empty N x N array
biphoton_spec = zeros((N,N), dtype = float)


for j in range(0, N):
    for i in range(0, N):

        # angular frequencies
        ws = ws_array[i]
        wi = wi_array[j]

        # refractive index
        ns = refractive_index_BBO(ws)
        ni = refractive_index_BBO(wi)
        np = refractive_index_BBO(wp)

        # Pump spectrum
        Ep_si = sqrt(Gaussian_spec(ws+wi, wp, Omega_p))

        # Alpha
        actual1 = alpha_func(ws, wi, ns, ni, Ep_si) #alpha_func2(ws+wi, ns+ni, Ep_si)#

        # Phi
        actual2 = Phi_L_func(L, wp, ws, wi, gamma, gamma, theta)

        # append to biphoton spectrum
        biphoton_spec[j,i] = actual1 * actual2

        if i == j:
            biphoton_spec[j,i] = 1.0

        if j == -i+N:
            biphoton_spec[j,i] = 1.0



#------------------------------------------------
#////////////////////////////////////////////////
#------------------------------------------------


#===============================================================================
#///////////////////////////////////////////////////////////////////////////////
#===============================================================================




#===============================================================================
# SIMULATION 2
#===============================================================================


#-------------------------
# Parameters (SI units)
#-------------------------
L       = 1e-2               # crystal length
theta   = -pi/4              # angle
Omega_p = 6e13               # proportional to the FWHM of the pump spectrum
gamma   = 8e-5 * 1e-12/1e-6  #

# pump central angular frequency spectrum
wp      = 2*pi*3e8/395e-9    # lambda_0 = 395e-9
#-----------------

#------------------------------------------------
# SIMULATION
#------------------------------------------------

# empty N x N array
biphoton_spec2 = zeros((N,N), dtype = float)


for j in range(0, N):
    for i in range(0, N):

        # angular frequencies
        ws = ws_array[j]
        wi = wi_array[i]

        # refractive index
        ns = refractive_index_BBO(ws)
        ni = refractive_index_BBO(wi)
        np = refractive_index_BBO(wp)

        # Ep
        Ep_si = sqrt(Gaussian_spec(ws+wi, wp, Omega_p))

        # alpha
        actual1 = alpha_func(ws, wi, ns, ni, Ep_si) #alpha_func2(ws+wi, ns+ni, Ep_si)#

        # Phi
        actual2 = Phi_L_func(L, wp, ws, wi, gamma, gamma, theta)

        # Append to biphoton spectrum
        biphoton_spec2[j,i] = actual1 * actual2

        if i == j:
            biphoton_spec2[j,i] = 1.0

        if j == -i+N:
            biphoton_spec2[j,i] = 1.0


#------------------------------------------------
#////////////////////////////////////////////////
#------------------------------------------------


#===============================================================================
#///////////////////////////////////////////////////////////////////////////////
#===============================================================================



#===============================================================================
# PLOTS
#===============================================================================


#---------------------------------------------------
# NEW COLORMAP
#---------------------------------------------------
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
bottom = cm.get_cmap('jet', 236)
top = cm.get_cmap('Blues', 20)
newcolors = vstack((top(linspace(0, 1, 20)),
                       bottom(linspace(0, 1, 236))))

newcmp = ListedColormap(newcolors, name='jet')
#---------------------------------------------------
#///////////////////////////////////////////////////
#---------------------------------------------------


#-----------------------------------
# PLOT ANGULAR FREQUENCY
#-----------------------------------
new_ws = ((ws_array-wp/2))/wp
new_wi = ((wi_array-wp/2))/wp

# figure
figure()
suptitle("Biphoton spectral amplitude "+r'$\left|A(\omega_s, \omega_i) \right|$')

# plot simulation 1
subplot(121)
title(r'$L=1\,\mathrm{cm},\,\,\,\,\theta=\pi/20$')
pcolormesh(new_wi, new_ws, abs(biphoton_spec), cmap=newcmp)
xlim(-0.04,0.04)
ylim(-0.04,0.04)
xlabel(r'$\left(\omega_s-\omega_p/2\right)/\omega_p$')
ylabel(r'$\left(\omega_i-\omega_p/2\right)/\omega_p$')
colorbar()

# plot simulation 2
subplot(122)
title(r'$L=1\,\mathrm{cm},\,\,\,\,\theta=-\pi/4$')
pcolormesh(new_wi, new_ws, abs(biphoton_spec2), cmap=newcmp)
xlim(-0.02,0.02)
ylim(-0.02,0.02)
xlabel(r'$\left(\omega_s-\omega_p/2\right)/\omega_p$')
ylabel(r'$\left(\omega_i-\omega_p/2\right)/\omega_p$')
colorbar()

#-----------------------------------
#//////////////////////////////////
#-----------------------------------



#-----------------------------------
# PLOT WAVELENGTH
#-----------------------------------

# wavelength array - signal
wave_array_s = 2*pi*3e8/ws_array
wave_array_s*=1e9
wave_array_s = flip(wave_array_s)


# wavelength array - idler
wave_array_i = 2*pi*3e8/wi_array
wave_array_i*=1e9
wave_array_i = flip(wave_array_i)

# PLOT
figure()
suptitle("Biphoton spectral amplitude "+r'$\left|A(\omega_s, \omega_i) \right|$')

# plot simulation 1
subplot(121)
title(r'$L=1\,\mathrm{cm},\,\,\,\,\theta=\pi/20$')
pcolormesh(wave_array_i, wave_array_s, abs(biphoton_spec), cmap=newcmp)
xlim(750,830)
ylim(750,830)
xlabel("signal wavelength (nm)")
ylabel("idler wavelength (nm)")
colorbar()

# plot simulation 2
subplot(122)
title(r'$L=1\,\mathrm{cm},\,\,\,\,\theta=-\pi/4$')
pcolormesh(wave_array_i, wave_array_s, abs(biphoton_spec2), cmap=newcmp)
xlim(750,830)
ylim(750,830)
xlabel("signal wavelength (nm)")
ylabel("idler wavelength (nm)")
colorbar()

#-----------------------------------
#//////////////////////////////////
#-----------------------------------


show()


#===============================================================================
#///////////////////////////////////////////////////////////////////////////////
#===============================================================================