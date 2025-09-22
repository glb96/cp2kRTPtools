#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Fourier Transform utils
# Contains functions used to perform time-fourier transform in 1D, for instance the total dipole moment, or in 3D, for instance the time-dependent electronic density.

import numpy as np
import os

au_to_ev = 27.211324570273 
au_to_Vm = 5.142206747E+11
au_to_D = (8.4783536255 / 3.335640952)
au_to_fs = 0.02418884254


#########################################################################
#########################################################################
#########################################################################

def creates_ft_authorized_omega(L_time, omega_min, omega_max, N_omega, warning_msg=True):
    '''
    We should perfrom the fourier transform in a discret nbr of point in time.
    To avoid problems, we have to ensure that: total_time*omega = N*2*pi
    Hence: omega can only be multiple of: 2*pi/total_time
    
    Note that the time and omega provided should have the same (inversed) unit!!!
    '''
    delta_omega = 2*np.pi/(L_time[-1]-L_time[0]) #in (a.u.)^(-1)
    omega_start = (omega_min // delta_omega) * delta_omega
    omega_end = ((omega_max // delta_omega) + 1) * delta_omega
    L_omega = np.arange(omega_start, omega_end + 0.5*delta_omega, delta_omega) # the frequency close to one desired which fullfies total_time*omega = N*2*pi
    if warning_msg:
        if len(L_omega) < N_omega:
            print('WARNING: you have asked ' + str(N_omega) + ' frequencies within:[' + str(omega_min) + ', ' + str(omega_max) + '] but the amount of time in this sample is too small to allow enough frequencies. Instead we can sample ' + str(len(L_omega)) + ' frequencies.')
        elif len(L_omega) > N_omega:
            print('WARNING: you have asked ' + str(N_omega) + ' frequencies within:[' + str(omega_min) + ', ' + str(omega_max) + '] but the amount of time in this sample is larger so that we can sample more frequencies. Instead we  sample ' + str(len(L_omega)) + ' frequencies.')
    return(L_omega)

#########################################################################
#########################################################################
#########################################################################

def fourier_transform(L_time, L_f, L_omega):
    '''
    Perform the fourier transform of something using the sin and cos definition.
    FFT can lead to some issue sometimes: better use this definition to better sleep.
    
    L_time and L_f should have the same length, and should be 1D.
    L_omega can have any length. 
    
    There are many ways of adding a smoothing of the time series. 
    Hence, if one want to perform a smoothing, it is easier to call this function with the 
    smoothed function directly. Note that some extra normalization may be needed if a damping is applied.
    '''
    N_omega = len(L_omega)
    N_t = len(L_time)
    delta_t = L_time[1]-L_time[0]
    L_ft = np.zeros(N_omega, dtype=complex)
    for k in range(0, N_omega, 1):
        omega = L_omega[k]
        for t in range(0, N_t, 1):
            L_ft[k] += L_f[t]*(np.cos(omega*L_time[t])+1j*np.sin(omega*L_time[t]))
    L_ft = L_ft/N_t # *delta_t/total_time = 1/N_t.
    return(L_ft)

#########################################################################
#########################################################################
#########################################################################

def exponential_damping(L_time, tau):
    '''
    Returns an exponential damp function according to the tau parameter
    '''
    return(np.exp((L_time[0]-L_time)/tau))

#########################################################################

def hann_damping(L_time, t0, sigma):
    '''
    Returns an hann damp function according to the t0 and sigma parameter.
    '''
    return(np.cos(np.pi/sigma*(L_time-t0))**2)

#########################################################################
#########################################################################
#########################################################################

def perform_ft(L_time, L_f, L_damp, omega_min, omega_max, N_omega, omega_precision=1):
    '''
    Perform Fourier Transform of a list L_f which can be of size NxM, which it associated time L_time of size N.
    omega_min and omega_max are the minimal and maximal frequencies, N_omega is the maximal number of frequencies probed.
    
    The time unit and the frequencies should be given in the same unit.
    
    Note that all the maximal number of frequencies required may not be used because of the conditions linking the frequencies and the total time.
    If you want to increase the number of frequencies probed within the desired range, you can use the omega_precision optional input. 
        + If set to 1 (default), the frequency are sample once
        
        + If set to 2, the total amount of available time will be slightly reduced in order to probe other frequencies in the required range. This can double the number of frequencies probed. 
        
        + If set to -1, the same procedure as described in the +2 option is repeated. The number of time this operation is repeated depends on the time step and the maximal frequency required: n_iter = int(2*np.pi/(delta_t*omega_max)). During this run, the total amount of time will be reduced by 1 time step, 2 time step.... n_iter time step. For each of these new total amount of time, the authorized frequencies are computed. The returned quantites are all the authorized frequencies and their associated complex amplitude. This procedure can be very slow, if it is the case I recommand to use omega_precision=2 instead or to reduce the [omega_min, omega_max] range.
    '''
    if omega_precision != 1 and omega_precision != 2 and omega_precision != -1:
         raise Exception('ERROR: omega_precision can only be equal to 1, 2 or -1.')
    
    shape = np.shape(L_f)
    if len(shape) == 1: #only 1D 
        M = 1
        L_f = np.array([L_f]).T
    elif len(shape) == 2:
        M = shape[1]
    else:
        raise Exception('ERROR: the shape of L_f must be MxN or N, not higher!')
    print('The dimension of the function to Fourier transform is:', M)
    
    if omega_precision == 1:  
        L_omega = creates_ft_authorized_omega(L_time, omega_min, omega_max, N_omega)
        L_ft = np.zeros((len(L_omega), M), dtype=complex)
    elif omega_precision == 2:    
        L_omega_1 = creates_ft_authorized_omega(L_time, omega_min, omega_max, N_omega)
        total_time = L_time[-1]-L_time[0]
        delta_t = L_time[1]-L_time[0]
        delta_omega = L_omega_1[1]-L_omega_1[0]
        shift_omega = L_omega_1[0]+delta_omega/2
        time_domain_shift = total_time%(2*np.pi/shift_omega)
        N_to_dump = int(-time_domain_shift//delta_t)
        L_time_2 = np.array(L_time[:N_to_dump])
        L_omega_2 = creates_ft_authorized_omega(L_time_2, omega_min, omega_max, N_omega)
        L_omega = np.append(L_omega_1, L_omega_2)
        ind =  np.argsort(L_omega, axis=0)
        L_omega = np.take_along_axis(L_omega, ind, axis=0)
        L_ft = np.zeros((len(L_omega), M), dtype=complex)
    elif omega_precision == -1:
        delta_t = L_time[1]-L_time[0]
        n_max_dump = int(2*np.pi/(delta_t*omega_max)) # the nbr of point to remove until we recover the same authorized frequency for the max frequency required
        print('The frequency range will be sampled ' + str(n_max_dump) + ' times.')
     
    
    for i in range(M):
        if omega_precision == 1:  
            L_ft[:,i] = fourier_transform(L_time,  L_f[:,i]*L_damp, L_omega)
        elif omega_precision == 2:  
            L_ft_1 = fourier_transform(L_time,  L_f[:,i]*L_damp, L_omega_1)
            L_ft_2 = fourier_transform(L_time_2,  L_f[:,i]*L_damp, L_omega_2)
            # merge the 2 FTs
            L_ft_i = np.append(L_ft_1, L_ft_2)
            L_ft[:, i] = np.take_along_axis(L_ft_i, ind, axis=0)
        elif omega_precision == -1:
            for n_dump in range(0, n_max_dump, 1):
                if n_dump == 0:
                    L_omega_temp = creates_ft_authorized_omega(L_time, omega_min, omega_max, N_omega, warning_msg=False)
                    L_ft_temp = fourier_transform(L_time,  L_f[:,i]*L_damp, L_omega_temp)
                    L_omega = np.array(L_omega_temp)
                    L_ft_i = np.array(L_ft_temp)
                else:
                    L_omega_temp = creates_ft_authorized_omega(L_time[:-n_dump], omega_min, omega_max, N_omega, warning_msg=False)
                    L_ft_temp = fourier_transform(L_time[:-n_dump],  L_f[:,i]*L_damp, L_omega_temp)
                    L_omega = np.append(L_omega, L_omega_temp)
                    ind =  np.argsort(L_omega, axis=0)
                    L_omega = np.take_along_axis(L_omega, ind, axis=0)
                    L_ft_i = np.append(L_ft_i, L_ft_temp)
                    L_ft_i = np.take_along_axis(L_ft_i, ind, axis=0)
            if i == 0:
                L_ft = np.zeros((len(L_omega), M), dtype=complex)
                L_ft[:, i]= np.array(L_ft_i)
            else:
                L_ft[:, i]= np.array(L_ft_i)
    return(L_omega, L_ft)

#########################################################################
#########################################################################
#########################################################################   

