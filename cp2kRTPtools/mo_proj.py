#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Projection utils:
# Perform the projection of a time-dependent state on reference one. 

import numpy as np
import os

#########################################################################
#########################################################################
#########################################################################

def compute_projection_ref_td(L_ref_MO, L_td_MO, L_ao_matrix=False):
    '''
    Project the Time-dependent MO into the reference one. 
    
    The reference MO can be from GS or a TD one, ie real or complex. 
    In any case, the projection is normalized with the amplitude of the reference state. 
    
    This function returns the amplitude and the phase of the projection.
    An amplitude of 1 means that the TD state and the reference one are the same. Note that this function does not use the overlap matrix: orthonormal states may have a non-zero projected value.
    The phase make sens if the projection is close to one: you can thus look at the evolution of the phase wrt time and extract the associated frequency of the state (its 'energy').
    '''
    nbr_ts =len(L_td_MO[:,0])
    L_projection = np.zeros(nbr_ts, dtype=complex)
    
    
    if isinstance(L_ao_matrix, bool):  # Approximate projection
        norm_ref = np.sum(np.absolute(L_ref_MO)**2)
        for TTT in range(0, nbr_ts, 1):
            L_projection[TTT] = np.sum(np.conjugate(L_ref_MO)*L_td_MO[TTT, :])
            
    else: # exact projection 
        nbr_ao = len(L_ao_matrix)
        norm_ref = 0
        for i in range(nbr_ao): 
            for j in range(nbr_ao):
                norm_ref += np.conjugate(L_ref_MO[i])*L_ao_matrix[i, j]*L_ref_MO[j]

        for TTT in range(0, nbr_ts, 1):
            for i in range(nbr_ao):
                for j in range(nbr_ao):
                    L_projection[TTT] += np.conjugate(L_ref_MO[i])*L_ao_matrix[i, j]*L_td_MO[TTT, j]
            
    L_projection = L_projection/norm_ref
        
    return(np.absolute(L_projection), np.angle(L_projection))

#########################################################################
#########################################################################
#########################################################################

