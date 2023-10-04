#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Open CP2K output files
# List of function used to open and read the output from CP2K software. 
# Please note that depending on the input parameters used, or the CP2K version, you may experience some issues while trying to open some file. So do not hesitate to update the function in this file: the rest of the python package should work rather independently. 

import numpy as np
import os

###########################################################################################
###########################################################################################
###########################################################################################

def read_one_line_energy(out):
    '''
    Read one line of the dipole moment file. 
    Return the dipole moment value. 
    Note that this function may have to be adapted if the cp2k output file format changes!
    '''
    out = out.replace('Total energy:',' ')
    out = float(out)
    return(out)

###########################################################################################

def read_energy(filename):
    temp_file_name = 'tempenerg.txt'
    temp_file_name2 = 'tempenerg2.txt'
    os.system('grep -B 4 "Time needed for propagation:" ' + filename + ' > ' + temp_file_name)
    os.system('grep "Total energy:" ' + temp_file_name + ' > ' + temp_file_name2)
    with open(temp_file_name2, 'r') as thefile:
        L_energy = []
        for line in thefile:
            L_energy.append(read_one_line_energy(line))
    L_energy = np.array(L_energy)
    os.system('rm ' + temp_file_name)
    os.system('rm ' + temp_file_name2)  
    return(L_energy)

###########################################################################################
###########################################################################################
###########################################################################################

def read_one_line_field(out):
    '''
    Read one line of the dipole moment file. 
    Return the dipole moment value. 
    Note that this function may have to be adapted if the cp2k output file format changes!
    '''
    out = out.replace('X=', '/')
    out = out.replace('Y=', '/')
    out = out.replace('Z=', '/')
    out = out.split('/')
    
    # This is needed because for very small field values, cp2k sometime write it weirdly: for instance like 0.47089303-100 instead of 0.47089303E-100 because it lacks character space. In the futur the best way would be to avoid writting such small value from CP2K directly (wrtting zero instead of 1E-100 ...).
    L_field_t = [0.0, 0.0, 0.0]
    for k in range(1, 4, 1):
        try:
            L_field_t[k-1] = float(out[k])
        except:
            pass
    return(L_field_t)

###########################################################################################

def read_applied_field(filename):
    temp_file_name = 'tempfield.txt'
    os.system('grep "X=" ' + filename + ' > ' + temp_file_name)
    with open(temp_file_name, 'r') as thefile:
        L_field = []
        for line in thefile:
            L_field.append(read_one_line_field(line))
    L_field = np.array(L_field)
    os.system('rm ' + temp_file_name)
    return(L_field)

###########################################################################################
###########################################################################################
###########################################################################################

def read_one_line_dipole_moment(out):
    '''
    Read one line of the dipole moment file. 
    Return the dipole moment value. 
    Note that this function may have to be adapted if the cp2k output file format changes!
    '''
    out = out.replace('X=', '/')
    out = out.replace('Y=', '/')
    out = out.replace('Z=', '/')
    out = out.replace('Total=', '/')
    out = out.split('/')
    L_dipole_t = [float(out[item]) for item in range(1, 4, 1)]
    return(L_dipole_t)

###########################################################################################

def read_dipole_moment(filename):
    temp_file_name = 'tempdipole.txt'
    temp_file_name2 = 'tempdipole2.txt'
    os.system('grep -A 1 "Dipole moment" ' + filename + ' > ' + temp_file_name)
    os.system('grep "X=" ' + temp_file_name + ' > ' + temp_file_name2)
    with open(temp_file_name2, 'r') as thefile:
        L_dipole = []
        for line in thefile:
            L_dipole.append(read_one_line_dipole_moment(line))
    L_dipole = np.array(L_dipole)
    os.system('rm ' + temp_file_name)
    os.system('rm ' + temp_file_name2)
    return(L_dipole)

###########################################################################################
###########################################################################################
###########################################################################################

def read_one_line_momentum(out):
    '''
    Read one line of the dipole moment file. 
    Return the dipole moment value. 
    Note that this function may have to be adapted if the cp2k output file format changes!
    '''
    out = out.replace('X=', '/')
    out = out.replace('Y=', '/')
    out = out.replace('Z=', '/')
    out = out.replace('Total=', '/')
    out = out.split('/')
    L_moment_t = [float(out[item]) for item in range(1, 4, 1)]
    return(L_moment_t)

###########################################################################################

def read_momentum(filename):
    temp_file_name = 'tempmomentum.txt'
    temp_file_name2 = 'tempmomentum2.txt'
    os.system('grep -A 1 "momentum operator" ' + filename + ' > ' + temp_file_name)
    os.system('grep "X=" ' + temp_file_name + ' > ' + temp_file_name2)
    with open(temp_file_name2, 'r') as thefile:
        L_moment = []
        for line in thefile:
            L_moment.append(read_one_line_momentum(line))
    L_moment = np.array(L_moment)
    os.system('rm ' + temp_file_name)
    os.system('rm ' + temp_file_name2)
    return(L_moment)

#########################################################################
#########################################################################
#########################################################################

def read_one_line_NL_terms_1(out):
    '''
    Read one line of the dipole moment file. 
    Return the dipole moment value. 
    Note that this function may have to be adapted if the cp2k output file format changes!
    '''
    out = out.replace('X=', '/')
    out = out.replace('Y=', '/')
    out = out.replace('Z=', '/')
    out = out.replace('Total=', '/')
    out = out.split('/')
    L_out_t = [float(out[item]) for item in range(1, 4, 1)]
    return(L_out_t)

###########################################################################################

def read_one_line_NL_terms_2(out):
    '''
    Read one line of the dipole moment file. 
    Return the dipole moment value. 
    Note that this function may have to be adapted if the cp2k output file format changes!
    '''
    out = out.replace('XX=', '/')
    out = out.replace('XY=', '/')
    out = out.replace('XZ=', '/')
    out = out.split('/')
    L_out_t = [float(out[item]) for item in range(1, 4, 1)]
    return(L_out_t)

###########################################################################################

def read_one_line_NL_terms_3(out):
    '''
    Read one line of the dipole moment file. 
    Return the dipole moment value. 
    Note that this function may have to be adapted if the cp2k output file format changes!
    '''
    out = out.replace('YY=', '/')
    out = out.replace('YZ=', '/')
    out = out.replace('ZZ=', '/')
    out = out.split('/')
    L_out_t = [float(out[item]) for item in range(1, 4, 1)]
    return(L_out_t)

###########################################################################################

def read_NL_terms(filename):
    '''
    Read one line of the dipole moment file. 
    Return the terms associated to non-local pseudo potential for the generalized current calculation in the velocity gauge.
    3 expectation values are read: 
    <[r, V_nl]>(t),  <r x V_nl x r>(t), and, <r x r x V_nl + V_nl x r x r>(t)
    
    Note that this function may have to be adapted if the cp2k output file format changes!
    '''
    # [r, V_nl]
    temp_file_name = 'temp.txt'
    temp_file_name2 = 'temp2.txt'
    os.system('grep -FA 1 "[r,V_nl]" ' + filename + ' > ' + temp_file_name)
    os.system('grep "X=" ' + temp_file_name + ' > ' + temp_file_name2)
    with open(temp_file_name2, 'r') as thefile:
        L_com_r_v = []
        for line in thefile:
            L_com_r_v.append(read_one_line_NL_terms_1(line))
    L_com_r_v = np.array(L_com_r_v)
    
    # r x V_nl x r
    temp_file_name = 'temp.txt'
    temp_file_name2 = 'temp2.txt'
    os.system('grep -FA 3 "r x V_nl x r" ' + filename + ' > ' + temp_file_name)
    os.system('grep "XX=" ' + temp_file_name + ' > ' + temp_file_name2)
    with open(temp_file_name2, 'r') as thefile:
        L_rvr = []
        for line in thefile:
            L_rvr.append(read_one_line_NL_terms_2(line))
    os.system('grep "YY=" ' + temp_file_name + ' > ' + temp_file_name2)
    with open(temp_file_name2, 'r') as thefile:
        n = 0
        for line in thefile:
            for element in read_one_line_NL_terms_3(line):
                L_rvr[n].append(element)
            n += 1
    L_rvr = np.array(L_rvr)
    
    # rx r x V_nl + V_nl x r x r
    temp_file_name = 'temp.txt'
    temp_file_name2 = 'temp2.txt'
    os.system('grep -FA 3 "r x r x V_nl + V_nl x r x r" ' + filename + ' > ' + temp_file_name)
    os.system('grep "XX=" ' + temp_file_name + ' > ' + temp_file_name2)
    with open(temp_file_name2, 'r') as thefile:
        L_rrv_vrr = []
        for line in thefile:
            L_rrv_vrr.append(read_one_line_NL_terms_2(line))
    os.system('grep "YY=" ' + temp_file_name + ' > ' + temp_file_name2)
    with open(temp_file_name2, 'r') as thefile:
        n = 0
        for line in thefile:
            for element in read_one_line_NL_terms_3(line):
                L_rrv_vrr[n].append(element)
            n += 1
    L_rrv_vrr = np.array(L_rrv_vrr)
    
    os.system('rm ' + temp_file_name)
    os.system('rm ' + temp_file_name2)
    
    return(L_com_r_v, L_rvr, L_rrv_vrr)

#########################################################################
#########################################################################
#########################################################################

def return_list_of_float_from_string_mo(mystring):
    '''
    Return from a list of number written in string the list of the float value in a list object. 
    If the last number added is a integer, return is_int=True, otherwise, is_int=False
    '''
    mylist = []
    to_add = False
    temp = ''
    is_int = False
    mystring = mystring.rstrip()
    for item in mystring:
        if item == ' ':
            if to_add:
                to_add = False
                mylist.append(float(temp))
                temp = ''
        else:
            to_add = True
            temp += item
    if len(temp):
        is_int = temp.isdigit()
        mylist.append(float(temp))
    return(mylist, is_int)

#########################################################################

def read_one_mo_set(L_one_mo_set):
    '''
    Return the AO coefficient for all the MO contains in the list of string. 
    
    In a standard CP2K MO output file, the MO are written 4 by 4 maximum. 
    You have the list of the AO coefficient for 4 MOs, then the list of the AO 
    coefficient for the next 4 MOs ect... 
    
    This function read one 'mo set' (4 MO maximum) at a time. 
    The output list size is nbr_MO X nbr_AO 
    '''
    L_mo = []
    L_new_mo = [[]]
    nbr_line = len(L_one_mo_set)
    shall_i_change_mo_subset = True
    shall_i_read_ao = False
    for ttt in range(1, nbr_line, 1):
        if shall_i_change_mo_subset:
            if len(L_new_mo[0]):
                for new_mo in L_new_mo:
                    L_mo.append(new_mo)
            L_list_of_mo_names, _ = return_list_of_float_from_string_mo(L_one_mo_set[ttt][4:])
            L_new_mo = [[] for k in range(len(L_list_of_mo_names))]
            shall_i_change_mo_subset = False
            shall_i_read_ao = False
        else:
            if L_one_mo_set[ttt][0:12] == ' MO|    1   ':
                shall_i_read_ao = True
            
            if shall_i_read_ao: 
                L_list_of_ao, is_int = return_list_of_float_from_string_mo(L_one_mo_set[ttt][29:])
                if len(L_list_of_ao):
                    if is_int: #new set of MOs
                        shall_i_change_mo_subset = True
                        shall_i_read_ao = False
                    else:
                        for k in range(len(L_list_of_mo_names)):
                            L_new_mo[k].append(L_list_of_ao[k])
    for new_mo in L_new_mo: 
        L_mo.append(new_mo)
    return(L_mo)

#########################################################################

def read_ref_mo(filename, spin_nbr=2):
    '''
    Read the MO coefficient at the end of an SCF static calculation. 
    Surprisingly, the name of the file containing the MO is given by the input 'filename'.
    By default, the number of spin is 2.
    You can use the optional input spin_nbr to adjust. 
    
    The size of the output is nbr_spin X nbr_MO X nbr_AO
    '''
    L_one_mo_set = []
    L_ref_mo = []
    spin_counter = 0
    shall_i_continue = True
    shall_i_read_mo_set = False
    with open(filename, 'r') as file:
        while shall_i_continue:
            L_value = file.readline()
            if 'EIGENVALUES' in L_value and not 'SCF' in L_value:
                shall_i_read_mo_set = True
                L_one_mo_set = []
            elif shall_i_read_mo_set:
                if not 'a.u.' in L_value and not 'eV' in L_value :
                    L_one_mo_set.append(L_value)
                else:
                    if len(L_one_mo_set) != 0:
                        L_mo = read_one_mo_set(L_one_mo_set)
                        L_ref_mo.append(L_mo)
                        L_one_mo_set = []
                        spin_counter += 1
                        if spin_counter == spin_nbr:
                            shall_i_continue = False
            else:
                pass
    return(np.array(L_ref_mo))
                
#########################################################################

def read_td_mo(filename, spin_nbr=4):
    '''
    Read the Time Dependent MO coefficient from a TD calculation.
    Surprisingly, the name of the file containing the MO is given by the input 'filename'.
    By default, the number of spin is 4: the real and imaginary part of the 2 possible spins. 
    
    You can use the optional input spin_nbr to adjust: use either 2 or 4. 
    
    The size of the output is nbr_time_step X nbr_MO X nbr_AO
    
    Note that the output is always for the alpha and beta spin.
    If only one spin, the beta spin list is empty.
    '''
    if spin_nbr != 2 and spin_nbr != 4:
        raise Exception('WARNING: spin number can only be 2 or 4. (there are both the immaginary and real part for each MO!!!')
        
    L_one_mo_set = []
    L_one_td_mo = []
    L_td_mo = []
    spin_counter = 0
    shall_i_continue = True
    shall_i_read_mo_set = False
    time_step = 1 
    stop_me = 0 
    with open(filename, 'r') as file:
        while stop_me < 10 and shall_i_continue:
            try:
                L_value = file.readline()
            except:
                shall_i_continue = False
            if not len(L_value):
                stop_me += 1
            elif shall_i_continue:
                stop_me = 0
                #if 'EIGENVALUES' in L_value and 'RTP' in L_value:
                if 'RTP' in L_value:
                    shall_i_read_mo_set = True
                    
                if 'RTP' in L_value and shall_i_read_mo_set:
                    if len(L_one_mo_set) != 0:
                        L_mo = read_one_mo_set(L_one_mo_set)
                        L_one_td_mo.append(L_mo)
                        spin_counter += 1
                        if spin_counter == spin_nbr:
                            end_set = True
                            L_td_mo.append(L_one_td_mo)
                            L_one_td_mo = []
                            spin_counter = 0
                            print('Read time step: ', time_step)
                            time_step += 1
                    L_one_mo_set = []
                elif shall_i_read_mo_set:
                    if not 'a.u.' in L_value and not 'eV' in L_value :
                        L_one_mo_set.append(L_value)
                else:
                    pass
            else:
                pass
    nbr_ts = len(L_td_mo[:])
    print('nbr_ts=', nbr_ts)
    n_ao = len(L_td_mo[0][0][0][:])
    print('n_ao=', n_ao)
    
    if spin_nbr == 2:
        n_mo_alpha = len(L_td_mo[0][0])
        print('n_mo=', n_mo_alpha)
        L_MO = np.zeros((nbr_ts, n_mo_alpha, n_ao), dtype=complex)
        for ttt in range(0, nbr_ts, 1):
            for ao in range(0, n_ao, 1):
                for mo in range(0, n_mo_alpha, 1):
                    L_MO[ttt, mo, ao] = L_td_mo[ttt][0][mo][ao] + 1j*L_td_mo[ttt][1][mo][ao]
        return(L_MO, [])
    elif spin_nbr == 4:
        n_mo_alpha = len(L_td_mo[0][0])
        n_mo_beta = len(L_td_mo[0][1])
        print('n_mo_alpha=', n_mo_alpha)
        print('n_mo_beta=', n_mo_beta)
        L_MO_alpha = np.zeros((nbr_ts, n_mo_alpha, n_ao), dtype=complex)
        L_MO_beta = np.zeros((nbr_ts, n_mo_beta, n_ao), dtype=complex)
        for ttt in range(0, nbr_ts, 1):
            for ao in range(0, n_ao, 1):
                for mo in range(0, n_mo_alpha, 1):
                    L_MO_alpha[ttt, mo, ao] = L_td_mo[ttt][0][mo][ao] + 1j*L_td_mo[ttt][1][mo][ao]
                for mo in range(0, n_mo_beta, 1):
                    L_MO_beta[ttt, mo, ao] = L_td_mo[ttt][2][mo][ao] + 1j*L_td_mo[ttt][3][mo][ao]
        return(L_MO_alpha, L_MO_beta)
    else:
        raise Exception('WARNING: code issue, should not happens')
    
#########################################################################
#########################################################################
#########################################################################

def return_list_of_float_from_string_ao(mystring):
    '''
    Return from a list of number written in string the list of the float value in a list object. 
    If the last number added is a integer, return is_int=True, otherwise, is_int=False
    '''
    
    mystring = mystring.rstrip('\n')
    
    if len(mystring) == 0:
        return([], False)
    
    mylist = []
    is_int = False
    mystring = mystring.split(' ')
    mystring = [item for item in mystring if item != '']
    for k in range(-1, -len(mystring)-1, -1):
        try:
            temp = float(mystring[k])
            mylist.append(temp)
        except:
            break
            
    if len(mylist):
        mylist.reverse()
        is_int = mystring[-1].isdigit()
    return(mylist, is_int)

#########################################################################

def fill_matrix_element(mylist, L_one_ao_set):
    '''
    Fill the set of AO elements with the read value.
    '''
    
    if len(mylist) == 0: # there are blank line sometimes
        pass
    else:
        if len(L_one_ao_set) == 0:
            for value in mylist:
                L_one_ao_set.append([value])
        else:
            for t in range(len(mylist)):
                L_one_ao_set[t].append(mylist[t])     
    
#########################################################################

def read_ao_matrix(filename, startkey):
    '''
    Read the AO matrix coefficient from filename.
    The function use the startkey string to start looking for the coefficient.
    
    The size of the output is nbr_AO X nbr_AO
    
    WARNING: note that this function may fail in the case where the ao_matrix is in a larger file and the matrix is small (ie less than 4x4).
    
    '''
    L_one_ao_set = []
    L_ao_matrix = []
    shall_i_read_mo_set = False
    line_read_from_file = 0 
    nbr_read_line = 0
    with open(filename, 'r') as file:
        total_nbr_line = sum([1 for line in file])
    
    with open(filename, 'r') as file:
        for t_read in range(0, total_nbr_line, 1):
            L_value = file.readline()
            if  startkey in L_value: # start recording after the startkey patern found
                shall_i_read_mo_set = True
                break 
                
        line_read_from_file = t_read

        if not shall_i_read_mo_set:
            raise Exception('ERROR: the function has not been able to find the patern ' + startkey + ' in the file ' + filename + ' in order to start reading the AO matrix.' )
        
        for t_read in range(line_read_from_file+1, total_nbr_line, 1): # read the matrix 
            L_value = file.readline()
            mylist, is_int = return_list_of_float_from_string_ao(L_value)
            if is_int: # start of a new lines 
                if len(L_ao_matrix) == 0 and len(L_one_ao_set) == 0: #  initialization
                    pass
                elif len(L_ao_matrix) == 0 and len(L_one_ao_set) != 0: # end of the initialization
                    L_one_ao_set = np.array(L_one_ao_set)
                    nbr_line, nbr_ao = L_one_ao_set.shape
                    L_ao_matrix = np.zeros((nbr_ao, nbr_ao))
                    for i in range(0, nbr_line, 1):
                        L_ao_matrix[nbr_read_line, :] = L_one_ao_set[i, :]
                        nbr_read_line += 1
                    L_one_ao_set = []    
                    
                    break   
                else:
                    raise Exception('CODE ERROR: this case should not happen!')
            else:
                fill_matrix_element(mylist, L_one_ao_set)
                
        line_read_from_file = t_read
        
        if len(L_one_ao_set) != 0: # the file has ended before the ao_matrix has been updated, this means that all the matrix fits on 4 column
            L_one_ao_set = np.array(L_one_ao_set)
            nbr_line, nbr_ao = L_one_ao_set.shape
            L_ao_matrix = np.zeros((nbr_ao, nbr_ao))
            for i in range(0, nbr_line, 1):
                L_ao_matrix[nbr_read_line, :] = L_one_ao_set[i, :]
                nbr_read_line += 1
            L_one_ao_set = []
            
        if nbr_read_line == nbr_ao: # all the lines have been read
            return(L_ao_matrix)
        
        for t_read in range(line_read_from_file+1, total_nbr_line, 1): # read the matrix 
            L_value = file.readline()
            mylist, is_int = return_list_of_float_from_string_ao(L_value)
            if is_int: # start of a new lines 
                pass          
            else: # fill several lines with column value
                fill_matrix_element(mylist, L_one_ao_set) 
                #print(len(L_one_ao_set[0]))
                if len(L_one_ao_set) and len(L_one_ao_set[0]) == nbr_ao:
                    L_one_ao_set = np.array(L_one_ao_set)
                    nbr_line, _ = L_one_ao_set.shape
                    for i in range(0, nbr_line, 1):
                        L_ao_matrix[nbr_read_line, :] = L_one_ao_set[i, :]
                        nbr_read_line += 1
                    L_one_ao_set = []
                    if nbr_read_line == nbr_ao:
                        break
                        
    if nbr_read_line == nbr_ao:
        return(L_ao_matrix)
    else:
        raise Exception('ERROR: the function was not able to read the whole AO matrix!')
        
#########################################################################
#########################################################################
#########################################################################

def read_mo_projection(filename):
    '''
    Read the population and the phase of projected Time-Dependent MO generated by FORCE_EVAL%DFT%REAL_TIME_PROPAGATION%PRINT%PROJECTION_MO.
    
    The output depends on the type of output requested during the cp2k run: the phase may not be present
    '''
    # Check the print output requested:
    temp_file_name = 'temp_projection.txt'
    case = 0 
    os.system('grep  "For each TD MOs required is printed: Population" ' + filename + ' > ' + temp_file_name)
    with open(temp_file_name, 'r') as thefile:
        thefile_read = [line for line in thefile]
        if len(thefile_read):
            nbr_time_step = len(thefile_read)
            if "Phase" in thefile_read[0]:
                case = 1 # several TD MO projection to read and phase 
            else:
                case = 2 # several TD MO projection to read and only the population
            
    os.system('grep  "The sum over all the TD MOs population:" ' + filename + ' > ' + temp_file_name)
    with open(temp_file_name, 'r') as thefile:
        thefile_read = [line for line in thefile]
        if len(thefile_read):
            nbr_time_step = len(thefile_read)
            case = 3 # several TD MO projection to read and only the population
    
    os.system('rm ' + temp_file_name)
    
    # call the right function depending on the type of output:
    if case == 1 or case == 2:
        return(read_mo_projection_several_td(filename, case, nbr_time_step))
    elif case == 3:
        return(read_mo_projection_sum_all_td(filename, nbr_time_step))
    else:
        raise Exception('ERROR: the type of output is not recognized of the file with filename ' + filename + '.')

#########################################################################

def read_mo_projection_several_td(filename, case, nbr_time_step):
    '''
    Read the projecton file in the cases where several TD MO are projected
    '''
    if case != 1 and case != 2:
        raise Exception('ERROR: this function should not be called in this case...')
    
    # determine how many MOs have to be read
    with open(filename, 'r') as thefile:
        thefile_read = [line for line in thefile]
        trotter = 0
        start_reading = False
        while trotter < len(thefile_read):
            myline = thefile_read[trotter].replace('\n', '')
            if "For each TD MOs required is printed: Population" in myline:
                if not start_reading: 
                    start_reading = True
                    n_td_mo = 0 
            else:
                if start_reading:
                    if len(myline) > 0 :
                        n_td_mo += 1
                    else:
                        start_reading = False
                        trotter = len(thefile_read) + 1
            trotter += 1
    
    # read and return
    if case == 1:
        L_popu = np.zeros((nbr_time_step, n_td_mo))
        L_phase = np.zeros((nbr_time_step, n_td_mo))
    elif case == 2:
        L_popu = np.zeros((nbr_time_step, n_td_mo))

    temp_file_name = 'temp_projection.txt'
    os.system('grep -A ' + str(n_td_mo) + ' "For each TD MOs required is printed: Population" ' + filename + ' > ' + temp_file_name)
    
    with open(filename, 'r') as thefile:
        trotter_t = 0
        trotter_mo = 0
        for line in thefile:
            if "For each TD MOs required is printed: Population" in line:
                start_reading = True
            else:
                if start_reading:
                    if trotter_mo == n_td_mo:
                        trotter_mo = 0
                        trotter_t += 1
                        start_reading = False
                    else:
                        if case == 1:
                            L_popu[trotter_t, trotter_mo], L_phase[trotter_t, trotter_mo] = read_mo_projection_popu(line, is_phase=True)
                        elif case == 2:
                             L_popu[trotter_t, trotter_mo] = read_mo_projection_popu(line)
                        trotter_mo += 1
                        
    os.system('rm ' + temp_file_name)
    
    if case == 1: 
        return(L_popu, L_phase)
    elif case == 2:
        return(L_popu)

#########################################################################

def read_mo_projection_sum_all_td(filename, nbr_time_step):
    '''
    Read the projecton file in the cases where several all the TD MO projection are summed up.
    '''
    
    temp_file_name = 'temp_projection.txt'
    os.system('grep -r "The sum over all the TD MOs population:" ' + filename + ' > ' + temp_file_name)
    L_popu = np.zeros(nbr_time_step)
    
    with open(temp_file_name, 'r') as thefile:
        thefile_read = [line for line in thefile]
        for t in range(nbr_time_step):
            myline = thefile_read[t].replace("The sum over all the TD MOs population:", "")
            myline = myline.replace("\n", "")
            L_popu[t] = float(myline)
    
    os.system("rm " + temp_file_name)
    
    return(L_popu)

#########################################################################

def read_mo_projection_popu(line, is_phase=False):
    '''
    Read the first and second float number of the line
    '''
    out = line.split(' ')
    elem = []
    for item in out:
        if len(item) > 0 :
            elem.append(float(item))
    if is_phase:
        return(elem[0], elem[1])
    else:
        return(elem[0])
    
#########################################################################
#########################################################################
#########################################################################



        









