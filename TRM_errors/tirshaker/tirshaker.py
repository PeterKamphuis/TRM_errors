# -*- coding: future_fstrings -*-

#Program to run tirshaker
from omegaconf import OmegaConf
import os
import copy
import psutil as psu
import random
import time
import re
from scipy import stats

from astropy.io import fits
import numpy as np
from dataclasses import dataclass, field
from typing import List
from TRM_errors.common.common import print_log,write_tirific,load_tirific,\
            set_format,set_limits,tirific_template,finish_current_run,check_pid
import subprocess
from datetime import datetime

class DefFileError(Exception):
    pass

class TirificOutputError(Exception):
    pass

class TirshakerInputError(Exception):
    pass

f'''
 NOTE:
 This is a re-imaged version of the code developed by G.I.G. Jozsa found at  https://github.com/gigjozsa/tirshaker.git

 Takes a tirific def file filename and varies it iterations times
 and runs it as many times to then calculate the mean and the
 standard deviation of parameters that have been varied, both are
 put into the tirific deffile outfilename.

 With parameter_groups parameter groups are defined that are varied
 homogeneously. This is a list of list of parameter names including
 the '='-sign.  Caution! At this stage it is essential that in the
 .def file there is the '='-sign directly attached to the parameter
 names and that there is a space between parameter values and the
 '='.

 For each parameter group (list of parameter names), the values of
 the rings specified in rings (which is a list of integers, the
 list member corresponding to the parameter with the same index in
 parameter_groups) are varied. Block is a list of indicators if the
 values should be varied by ring (similar to the !-sign in the VARY
 parameter of tirific) or as a whole, again indices indicating the
 associated parameter group. Variation quantifies the maximum
 variation of the parameter (again indices matching), with
 variation type (indices matching) 'a' indicating an absolute
 variation, 'r' indicating a relative one. Parameters are varied by
 a uniform variation with maximal amplitude variation. So if a
 parameter is x and v the variation, the parameter gets changed by
 a number between -v and v in case the matching variation type is
 'a' and it gets changed by a number between -v*x and v*x if the
 matching variation type is 'r'. Tirific is started iterations
 times with varied input parameters but otherwise identical
 parameters (except for output files which are suppressed) and the
 results are recorded. For each varied parameter the mean and
 standard deviation is calculated and returned in the output .def
 file outfilename. In outfilename LOOPS is set to 0 and any output
 parameter is preceded by a prefix outfileprefix. Random_seed is
 the random seed to make the process deterministic. If mode ==
 'mad', the median and the MAD is calculated and, based on that,
 values beyond a 3-sigma (sigma estimated from MAD) bracket around
 the median are rejected before calculating mean and standard
 deviation.

'''



'''Extract the existing errors if present'''
def get_existing_errors(Tirific_Template, fit_groups, log =False, verbose=True):
    log_statement = ''
    for group in fit_groups:
        fit_groups[group]['ERRORS'] = []
        for disk in fit_groups[group]['DISKS']:
            if disk == 1:
                par = group
            else:
                par = f"{group}_{disk}"
            errors = load_tirific(Tirific_Template,[f'# {par}_ERR'],array=True )
            if len(errors) == 0.:
                fit_groups[group]['ERRORS'].append(None)
            else:
                rings = fit_groups[group]['RINGS'][f'{disk}']
                fit_groups[group]['ERRORS'].append(np.mean(errors[rings[0]:rings[1]+1]))

    return log_statement
'''Extract the fitting parameters from PARMAX, PARMIN and DELSTART  '''
def get_variations(Tirific_Template, fit_groups, log =False, verbose =True):
    log_statement = ''
    parmax = Tirific_Template['PARMAX'].split()
    parmin = Tirific_Template['PARMIN'].split()
    delta_start = Tirific_Template['DELSTART'].split() 
    delta_end = Tirific_Template['DELEND'].split()
    for group in fit_groups:
        fit_groups[group]['PARMAX'] = float(parmax[fit_groups[group]['COLUMN_ID']])
        fit_groups[group]['PARMIN'] = float(parmin[fit_groups[group]['COLUMN_ID']])
        fit_groups[group]['FIT_DELTA'] = (float(delta_start[fit_groups[group]['COLUMN_ID']])+float(delta_end[fit_groups[group]['COLUMN_ID']]))/2.
    return log_statement

'''Extract the fitting parameters from the fitting VARY line '''
def get_groups(in_groups, var_index = {}, no_rings = 3, log = False, verbose=True):
    group_dict = {}   
    log_statement = ''
    if verbose:
        log_statement += print_log(f'''GET_GROUPS: We have found the following unformatted groups from VARY:
{'':8s}{in_groups}
''',log)
    for i,group in enumerate(in_groups):
        if verbose:
            log_statement += print_log(f'''GET_GROUPS: We are processing {group}
''',log)
        parameter = group.split()
        #first replace i with ! if i is present
        if  parameter[0][0] == 'i':
            parameter[0] = f'!{parameter[0][1:]}'
        #Tirific accepts spaces after the !
        if  parameter[0][0] == '!' and len(parameter[0]) == 1:
            parameter[0] = f'!{parameter[1]}'
            parameter.pop(1)
      
        count = 1
        #current_parameter = f'{re.sub("[^a-zA-Z]+", "", parameter[0])}_{count}'
        base_parameter =  f'{parameter[0].split("_")[0]}'
        if  base_parameter[0][0] == '!':
            base_parameter = base_parameter[1:]
        current_parameter = f'{base_parameter}_{count}'
        while current_parameter in group_dict:
            count += 1
            current_parameter =  f'{base_parameter}_{count}'
       
        group_dict[current_parameter] = {'COLUMN_ID': i}
        disks = parameter[0].split('_')
        try:
            group_dict[current_parameter]['DISKS'] =  [int(float(disks[1]))]  
            current_disk = disks[1]
        except IndexError:
            group_dict[current_parameter]['DISKS'] =  [1] 
            current_disk = '1'
        #Individual or block
        if parameter[0][0] == '!':
            group_dict[current_parameter]['BLOCK'] = False      
        else:
            group_dict[current_parameter]['BLOCK'] = True
        group_dict[current_parameter]['RINGS'] = {'EXTEND': []}

        # for singular parameters we do not have any rings
        if base_parameter in ['CONDISP']: 
            group_dict[current_parameter]['RINGS']['EXTEND'] =np.array([1,1])
            group_dict[current_parameter]['RINGS']['1'] = np.array([1,1])
        else:    
            start_ring = 0
            for i,part in enumerate(parameter[1:]):
                if part[0].isnumeric():
                    if start_ring == 0:
                        if ':' in part:
                            in_rings = [int(float(x)) for x in part.split(':')]
                            if in_rings.sort():
                                in_rings.sort()
                            start_ring = in_rings[0]
                        else:
                            start_ring = int(float(part))
                    last = False
                    try: 
                        if parameter[i+2][0].isnumeric():
                        #if the next batch is also numeri we continue to the end 
                            pass
                        else:
                            last =True
                    except IndexError:
                        last=True

                    if last:
                        if ':' in part:
                            in_rings = [int(float(x)) for x in part.split(':')]
                            if in_rings.sort():
                                in_rings.sort()
                            in_rings = np.array([start_ring,in_rings[-1]],dtype=int)      
                        else:
                            in_rings = np.array([start_ring,int(float(part))],dtype=int)

                        if current_disk not in group_dict[current_parameter]['RINGS']:
                            group_dict[current_parameter]['RINGS'][current_disk] = in_rings
                        else:
                            if verbose:
                                print(f'processing this group {group} for {current_parameter} and we have {in_rings} where we already have {group_dict[current_parameter]["RINGS"][current_disk]}')
                            raise DefFileError("The VARY settings in this deffile are not acceptable you have multiple indication of the same disk in one block.")
                        if len(group_dict[current_parameter]['RINGS']['EXTEND']) == 0:
                            group_dict[current_parameter]['RINGS']['EXTEND'] = in_rings
                        else:
                            group_dict[current_parameter]['RINGS']['EXTEND'][0] = \
                                np.nanmin([group_dict[current_parameter]['RINGS']['EXTEND'][0],in_rings[0]]) 
                            group_dict[current_parameter]['RINGS']['EXTEND'][1] = \
                                np.nanmax([group_dict[current_parameter]['RINGS']['EXTEND'][1],in_rings[1]])
                else:
                    start_ring = 0
                    disks = part.split('_')
                    try:
                        group_dict[current_parameter]['DISKS'].append(int(float(disks[1])))
                        current_disk = disks[1]
                    except IndexError:
                        group_dict[current_parameter]['DISKS'].append(1)  
                        current_disk = '1'
            if len(group_dict[current_parameter]['RINGS']['EXTEND']) == 0:
                group_dict[current_parameter]['RINGS']['EXTEND'] = np.array([1,no_rings],dtype=int)
            for disk in  group_dict[current_parameter]['DISKS']:           
                if f'{disk}' not in group_dict[current_parameter]['RINGS']:
                    group_dict[current_parameter]['RINGS'][f'{disk}'] = np.array([1,no_rings],dtype=int)

                if group_dict[current_parameter]['RINGS'][f'{disk}'][0] == group_dict[current_parameter]['RINGS'][f'{disk}'][1] \
                    and len(group_dict[current_parameter]['DISKS']) == 1:
                    group_dict[current_parameter]['BLOCK'] = False
        if verbose:  
            log_statement += print_log(f'''GET_FIT_GROUPS: We determined the group {group_dict[current_parameter]}
''',log) 
  
    return group_dict, log_statement

def set_fitted_variations(fit_groups,log=False,verbose=True):
    log_statement = ''
    for group in fit_groups:
        fit_groups[group]['VARIATION'] = [0., 'a']
        for i,disk in enumerate(fit_groups[group]['DISKS']):
            if fit_groups[group]['ERRORS'][i] == None:
                ini_var = 0.
            else: 
                ini_var = fit_groups[group]['ERRORS'][i]
            if fit_groups[group]['FIT_DELTA']*3. > ini_var:
                ini_var =  fit_groups[group]['FIT_DELTA']*3.
            if ini_var >  fit_groups[group]['VARIATION'][0]:
                fit_groups[group]['VARIATION'][0] = ini_var
        #fit_groups[group]['VARIATION'][0] *= 3.
    return log_statement

def set_manual_variations(fit_groups,variation= None,\
                                    cube_name= None,log=False,\
                                    verbose=True):
    log_statement = ''
    hdr = fits.getheader(cube_name)
    if not 'CUNIT3' in hdr:
        if abs(hdr['CDELT3']) > 100:
            hdr['CUNIT3'] = 'm/s'
        else:
            hdr['CUNIT3'] = 'km/s'
        if verbose:
            log_statement += print_log(f'''CLEAN_HEADER:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Your header did not have a unit for the third axis, that is bad policy.
{"":8s} We have set it to {hdr['CUNIT3']}. Please ensure that this is correct.'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
''',log)
    if hdr['CUNIT3'].upper() == 'HZ' or hdr['CTYPE3'].upper() == 'FREQ':
        if verbose:
            log_statement += print_log('CLEAN_HEADER: FREQUENCY IS NOT A SUPPORTED VELOCITY AXIS.', log)
        raise TirshakerInputError('The Cube has frequency as a velocity axis this is not supported')

    
    if hdr['CUNIT3'] == 'm/s':
        hdr['CDELT3'] = hdr['CDELT3']/1000.
   
    for group in fit_groups:
        if verbose:
            log_statement += print_log(f'''SET_MANUAL_VARIATIONS: processing {group}
''',log)
        groupbare = group.split('_')
        var_input = copy.deepcopy(getattr(variation,groupbare[0]))
        if var_input[1].lower() == 'res':
            if var_input[2].lower() == 'arcsec':
                var_input[0] = var_input[0] * hdr['BMAJ']*3600.
            elif var_input[2].lower() == 'degrees':
                var_input[0] = var_input[0] * hdr['BMAJ']
            elif var_input[2].lower() == 'angle':
                raise TirshakerInputError(f''' We have no way to relate an angle to the resolution of the cube''')
            elif var_input[2].lower() == 'km/s':
                var_input[0] = var_input[0] * hdr['CDELT3']
            elif var_input[2].lower() == 'm/s':
                var_input[0] = var_input[0] * hdr['CDELT3']*1000.
            elif var_input[2].lower() == 'Jy*km/s/arcsec^2':
                if 'NOISE' in hdr or 'FATNOISE' in hdr:
                    try:
                        noise = hdr['FATNOISE']
                    except KeyError:
                        noise = hdr['NOISE']
                else:
                   raise TirshakerInputError(f''' We have no way to relate an SBR to the cube without the noise level in the header''') 
                noise = noise*(2. *np.pi / (np.log(256.)))\
                    *hdr['BMAJ']*hdr['BMIN']*3600**2
                var_input[0] = var_input[0] * hdr['CDELT3']*noise
        fit_groups[group]['VARIATION'] = [var_input[0], var_input[3]]
    return log_statement


def get_manual_groups(cfg, rings = 1, cube_name='None', log= False, verbose= True):
    #First we get the groups that were fitted from file
    groups = cfg.variations.VARY.split(',')
    # And lets translate to a dictionary with the various fitting parameter type and
    # first we disentangle the tirific fitting syntax into a dictionary
    fit_groups,log_statement = get_groups(groups, no_rings = rings, log = log,verbose=cfg.general.verbose)
    #Then we set the  variation we want for the tirshaker for every group

    log_statement += set_manual_variations(fit_groups,variation=cfg.variations,\
                                    cube_name=cube_name,log=log,verbose=cfg.general.verbose)
    return fit_groups


def get_fitted_groups(Tirific_Template, log= False,verbose=True):
    #First we get the groups that were fitted from file
    groups = Tirific_Template['VARY'].split(',')
    # And lets translate to a dictionary with the various fitting parameter type and
    # first we disentangle the tirific fitting syntax into a dictionary
    fit_groups,log_statement = get_groups(groups, log = log,verbose=verbose)
    # !!!! There is no need to check the varindex as  Tirific ignores those values and interpolates any way 
    # Then we attach the fiiting variations to the groups
    log_statement += get_variations(Tirific_Template, fit_groups, log=log,verbose=verbose)
    # Check whether there are any errors present in the def file
    log_statement += get_existing_errors(Tirific_Template, fit_groups,log=log,verbose=verbose)
    #Then we set the  variation we want for the tirshaker for every group
    log_statement += set_fitted_variations(fit_groups,log=log,verbose=verbose)
    return fit_groups

get_fitted_groups.__doc__ =f'''
 NAME:
    get_fit_groups
 PURPOSE:
    get the groups that are fitting, whether they are a block or not and their expected errors

 CATEGORY:
    support_functions

 INPUTS:
    Tirific_Template =  the def file to get errors.
 OPTIONAL INPUTS:

 KEYWORD PARAMETERS:

 OUTPUTS:
     tirshaker settings
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 EXAMPLE:

 NOTE:
'''

    
def run_tirific(current_run, work_dir = os.getcwd(),deffile = 'tirific.def' ,tirific_call= 'tirific',\
                log=False, max_ini_time = 600, verbose = True):
    log_statement = print_log(f'''Starting a tirific run
''',log)
    restart_file = f"{work_dir}/restart_Error_Shaker.txt"   
    #Get the output fits file and def file defined in workdir+ deffile
    output_deffile = load_tirific(f'{work_dir}/{deffile}', Variables = ['TIRDEF'])[0]
    # Then if already running change restart file
    restart = False
    try:
        if check_pid(current_run.pid):
            restart = True
    except:
        pass 
    if restart:
        log_statement += print_log(f'''RUN_TIRIFIC: We are using an initialized tirific in {work_dir} with the file {deffile}
''',log)

        with open(restart_file,'a') as file:
            file.write("Restarting from previous run \n")
    else:
        log_statement += print_log(f'''RUN_TIRIFIC: We are initializing a new TiRiFiC in {work_dir} with the file {deffile}
''',log)
        with open(restart_file,'w') as file:
            file.write("Initialized a new run \n")
        current_run = subprocess.Popen([tirific_call,f"DEFFILE={deffile}","ACTION= 1"],\
                                       stdout = subprocess.PIPE, stderr = subprocess.PIPE,\
                                       cwd=work_dir,universal_newlines = True)

    currentloop =1
    max_loop = 0
    counter = 0

    current_process= psu.Process(current_run.pid)

    initialized = datetime.now()
    '''
    if Configuration['TIMING']:
        time.sleep(0.1)
        with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt",'a') as file:
            file.write(f"# TIRIFIC: Initializing Tirific at stage = {fit_type}, Loop = {Configuration['ITERATIONS']} {datetime.now()} \n")
            CPU,mem = get_usage_statistics(Configuration,current_process)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb for TiRiFiC \n")
    else:
    '''
    time.sleep(0.1)

    if verbose:
        print(f"\r{'':8s}RUN_TIRIFIC: 0 % Completed", end =" ",flush = True)
    triggered = False
    for tir_out_line in current_run.stdout:
        tmp = re.split(r"[/: ]+",tir_out_line.strip())
        counter += 1
        '''
        if (counter % 50) == 0:
            if Configuration['TIMING']:
                with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt",'a') as file:
                    if tmp[0] == 'L' and not triggered:
                        if tmp[1] == '1':
                            file.write(f"# TIRIFIC: Started the actual fitting {datetime.now()} \n")
                    CPU,mem = get_usage_statistics(Configuration,current_process)
                    file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb for TiRiFiC \n")
        '''
        if tmp[0] == 'L':
            if not triggered:
                triggered = True
            if int(float(tmp[1])) != currentloop and verbose:
                print(f"\r{'':8s}RUN_TIRIFIC: {set_limits(float(tmp[1])-1.,0.,float(max_loop))/float(max_loop)*100.:.1f} % Completed", end =" ",flush = True)
            currentloop  = int(float(tmp[1]))
            if max_loop == 0:
                max_loop = int(float(tmp[2]))
 

        if tmp[0].strip() == 'Finished':
            break
        if tmp[0].strip() == 'Abort':
            break
        if not triggered:
            #Check that the initialization doesn't take to long
            check = datetime.now()
            diff = (check-initialized).total_seconds()
            if diff > max_ini_time:
                log_statement += print_log(f'''RUN_TIRIFIC: After {diff/60.} min we could not find the expected output from the tirific run. 
running in the directory = {work_dir} 
and the file deffile = {deffile}                         
''',log)
                raise TirificOutputError(f'''{diff/60.} minutes after initialization the fitting has still not started.
We were running {deffile} and failed to find the output {output_deffile}.
''')
    if verbose:
        print(f'\n')
    '''
    if Configuration['TIMING']:
        with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt",'a') as file:
            file.write(f"# TIRIFIC: Finished this run {datetime.now()} \n")
            CPU,mem = get_usage_statistics(Configuration,current_process)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb for TiRiFiC \n")
    '''
    if verbose:
        print(f"{'':8s}RUN_TIRIFIC: Finished the current tirific run.")

    #The break off goes faster sometimes than the writing of the file so let's make sure it is present
    time.sleep(1.0)
    wait_counter = 0
    
    while not os.path.exists(f'{work_dir}/{output_deffile}') and wait_counter < 100.:
        time.sleep(0.3)
        wait_counter += 1
        if wait_counter/10. == int(wait_counter/10.):
            print(f"\r Waiting for {output_deffile}. \n", end = "", flush = True)
            log_statement += print_log(f'''RUN_TIRIFIC: we have waited {0.3*wait_counter} seconds for the output of tirific but it is not there yet.
''',log)
        if not  os.path.exists(output_deffile):
            log_statement += print_log(f'''RUN_TIRIFIC: After 30 seconds we could not find the expected output from the tirific run. We are raising an error for this galaxy.
''',log)
            raise TirificOutputError(f'''The tirific subprocess did not produce the correct output, most likely it crashed.
We were running {deffile} and failed to find the output  {output_deffile}.
''')
    if not log:
        return current_run
    else:
        return current_run,log_statement

run_tirific.__doc__ =f'''
 NAME:
    run_tirific

 PURPOSE:
    Check whether we have an initialized tirific if not initialize and run else restart the initialized run.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    current_run = subprocess structure of tirific

 OPTIONAL INPUTS:


    fit_type = 'Undefined'
    type of fitting

    stage = 'initial'
    stage of the fitting process

    max_ini_time = 600
    maximum time it can take for tirific to initialize 
    Higher ini times take longer

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_individual_iteration(Tirific_Template, i,fit_groups, directory,tirific_call, log = True,
                            name_in = 'Error_Shaker_In.def',name_out = 'Error_Shaker_Out.def', verbose=True,\
                            clean = True):
    log_statement = ''
    Current_Template = copy.deepcopy(Tirific_Template)
    Current_Template['RESTARTID']= i
    Current_Template['TIRDEF'] = name_out
    nur = int(float(Current_Template['NUR']))
    # Provide some info where we are
    if verbose:
        log_statement += print_log(f'''
    ******************************
    ******************************
    *** Tirshaker iteration {i:02d} ***
    ******************************
    ******************************
''',log)
    #Looping through all block
    for group in fit_groups:
        if group not in ['COLLECTED','TO_COLLECT']:
            if fit_groups[group]['BLOCK']:
                #If a block use the same variation for all rings in the groups
                variations = [fit_groups[group]['VARIATION'][0]*random.uniform(-1.,1.)]\
                    *(fit_groups[group]['RINGS']['EXTEND'][1]-fit_groups[group]['RINGS']['EXTEND'][0]+1)
            else:
                #If not a block use a different variation for all rings in the groups
                variations = [fit_groups[group]['VARIATION'][0]*random.uniform(-1.,1.) for x \
                            in range(fit_groups[group]['RINGS']['EXTEND'][0],fit_groups[group]['RINGS']['EXTEND'][1]+1)]

            for disk in fit_groups[group]['DISKS']:
                para = group.split('_')[0]
                if disk != 1:
                    para = f'{para}_{disk}'
                current_list = [float(x) for x in Current_Template[para].split()]
                while len(current_list) < nur:
                    current_list.append(current_list[-1])
                
                for l in range(fit_groups[group]['RINGS'][f'{disk}'][0],fit_groups[group]['RINGS'][f'{disk}'][1]+1):
                    if fit_groups[group]['VARIATION'][1] == 'a':
                        current_list[int(l-1)] += variations[int(l-fit_groups[group]['RINGS'][f'{disk}'][0])]
                    else:
                        current_list[int(l-1)] *= (1+variations[int(l-fit_groups[group]['RINGS'][f'{disk}'][0])])
                    #tirific does weird things when the initial values are outside parmax- pamin so check
                    if current_list[int(l-1)] <  fit_groups[group]['PARMIN']:
                        current_list[int(l-1)] = fit_groups[group]['PARMIN']+abs(current_list[int(l-1)]*0.05)
                    if current_list[int(l-1)] >  fit_groups[group]['PARMAX']:
                        current_list[int(l-1)] = fit_groups[group]['PARMAX']-abs(current_list[int(l-1)]*0.05)

                format = set_format(para)
                Current_Template[para] = ' '.join([f'{x:{format}}' for x in current_list])
    write_tirific(Current_Template, name =f'{directory}/{name_in}',full_name= True )
    output = {'i': i, 'directory': directory, 'deffile': name_in,\
             'tirific_call': tirific_call, 'TO_COLLECT':fit_groups['TO_COLLECT'], 'log': log_statement,\
             'tmp_name_out': Current_Template['TIRDEF'], 'verbose': verbose, 'clean': clean }
    return output

def run_individual_iteration(dict_input, log = True):
    current_run = None
    current_run = run_tirific(current_run,deffile=dict_input['deffile'],work_dir = dict_input['directory'],tirific_call=dict_input['tirific_call'], \
                            max_ini_time= 600, verbose=dict_input['verbose'])
    dict_input['log'] += finish_current_run(current_run,log=log)
        # Read the values of the pararameter groups
    output = {'log':dict_input['log']}    
    for parameter in dict_input['TO_COLLECT']:
            output[parameter] = load_tirific(f"{dict_input['directory']}/{dict_input['tmp_name_out']}",\
                    Variables = [parameter],array=True)
    # Remove the files after reading.    
    if dict_input['clean']:
        if os.path.isfile(f"{dict_input['directory']}/{dict_input['tmp_name_out']}"):
            os.remove(f"{dict_input['directory']}/{dict_input['tmp_name_out']}")
        if os.path.isfile(f"{dict_input['directory']}/{dict_input['deffile']}"):
            os.remove(f"{dict_input['directory']}/{dict_input['deffile']}")

    return output

def tirshaker_cleanup(fit_groups,cfg):
    log_statement = ''
    #read the original input
    Tirific_Template = tirific_template(filename=f'{cfg.general.directory}/{cfg.tirshaker.deffile_in}')

    fit_groups['FINAL_ERR'] = {}     
    for parameter in fit_groups['TO_COLLECT']:
        base_parameter= parameter.split('_')[0] 
        if cfg.general.verbose:
            log_statement += print_log(f'Processing {parameter}')
        all_iterations = np.array(fit_groups['COLLECTED'][parameter],dtype=float)
        original = load_tirific(Tirific_Template,\
                    Variables = [parameter],array=True)
        fit_groups['FINAL_ERR'][parameter] = np.zeros(all_iterations[0].size) 
        minimum_err = getattr(cfg.min_errors,base_parameter)
        for ring in range(all_iterations[0].size):
            all_its = all_iterations[:,ring]
            
            if cfg.general.calc_mode == 'mad':
                median = np.median(all_its)
                mad = stats.median_abs_deviation(all_its)
                madsigma = stats.median_abs_deviation(all_its) 
                average = np.average(all_its) 
                # Wow, np.std is the standard deviation using N and not N-1 in the denominator. So one has to use
                #std = np.sqrt(float(len(allparamsturned[j][k][l]))/float(len(allparamsturned[j][k][l])-1))*np.std(np.array(allparamsturned[j][k][l]))  
                std = np.std(all_its,ddof=1)     
                final = stats.tmean(all_its, (median-3*madsigma, median+3*madsigma))
                final_err =  stats.tstd(all_its, (median-3*madsigma, median+3*madsigma))
                if final_err > minimum_err:
                    fit_groups['FINAL_ERR'][parameter][ring] = final_err
                else:
                    fit_groups['FINAL_ERR'][parameter][ring] = minimum_err
                if cfg.general.verbose:
                    log_statement += print_log(f'TIRSHAKER: Parameter: {parameter} Ring: {ring} Pure average+-std: {average:.3e}+-{std:.3e} Median+-madsigma: {median:.3e}+-{madsigma:.3e} Average+-sigma filtered: {final:.3e}+-{final_err:.3e} \n')
            elif cfg.general.calc_mode == 'fat':
                #Fat does a lot of corrections which are not necessarily accounted for  
                #means that we should accont for the difference between the mean and the FAT value
                median = np.median(all_its)
                madsigma = stats.median_abs_deviation(all_its)
                mad_final = stats.tmean(all_its, (median-3*madsigma, median+3*madsigma))
                mad_final_err =  stats.tstd(all_its, (median-3*madsigma, median+3*madsigma))
                final = float((original[ring]+median)/2.)
                final_err = abs(np.sqrt((mad_final_err)**2+(abs(original[ring]-median)/2.)**2))
                if final_err > minimum_err:
                    fit_groups['FINAL_ERR'][parameter][ring] = final_err
                else:
                    fit_groups['FINAL_ERR'][parameter][ring] = minimum_err
                if cfg.general.verbose:
                    log_statement += print_log(f'TIRSHAKER: Parameter: {parameter} Ring: {ring} Final+-std: {final:.3e}+-{final_err:.3e} Median+-madsigma: {median:.3e}+-{madsigma:.3e} Average+-sigma filtered: {mad_final:.3e}+-{mad_final_err:.3e} \n')
 
          


    for parameter in fit_groups['TO_COLLECT']:
        format = set_format(parameter)
        Tirific_Template.insert(f'{parameter}',f'# {parameter}_ERR',f"{' '.join([f'{x:{format}}' for x in fit_groups['FINAL_ERR'][parameter]])}")     
    # Put them into the output file
    # Write it to a copy of the file replacing the parameters


   
    write_tirific(Tirific_Template, name = f'{cfg.general.directory}/{cfg.tirshaker.directory}/{cfg.tirshaker.deffile_out}',full_name=True)
    if cfg.general.verbose:
        log_statement += print_log(f'This is the final File with the errors {cfg.general.directory}/{cfg.tirshaker.directory}/{cfg.tirshaker.deffile_out}')
    return log_statement



def prepare_template(cfg, log=False,verbose=True):
    log_statement = ''
    #Read in the deffile
    Tirific_Template = tirific_template(filename=f'{cfg.general.directory}/{cfg.tirshaker.deffile_in}')
    # First we make a directory to keep all contained
    if not os.path.isdir(f'{cfg.general.directory}/{cfg.tirshaker.directory}/'):
        os.mkdir(f'{cfg.general.directory}/{cfg.tirshaker.directory}/')
    Tirific_Template['RESTARTNAME']= f"restart_Error_Shaker.txt"
   
    if cfg.tirshaker.individual_loops == -1:
        pass
    else:
        Tirific_Template['LOOPS'] = cfg.tirshaker.individual_loops
     #Some parameters are stripped by tirific so to make sure they are always there we add theif not present
    if 'GR_CONT' not in Tirific_Template:
        Tirific_Template['GR_CONT']=' '
    if 'RESTARTID' not in Tirific_Template:
        Tirific_Template['RESTARTID'] = '0'
    
    #Determine the error block from the  fit settings.
    if cfg.tirshaker.mode == 'fitted':
        fit_groups = get_fitted_groups(Tirific_Template,verbose=cfg.general.verbose)
    elif cfg.tirshaker.mode == 'manual':
        fit_groups = get_manual_groups(cfg, rings = int(float(Tirific_Template['NUR'])),\
                            cube_name = Tirific_Template['INSET'],verbose=cfg.general.verbose)
        #!!!!!!!!!!!!!!!!!!!!!!!!For this to wrok properly we should set the fitting parameters to the template
    else:
        if verbose:
            log_statement += print_log(f'''RUN_TIRSHAKER: The Tirshaker mode {cfg.tirshaker.mode} is not yet fully functional. Please use a different mode
''',log)
        raise TirshakerInputError(f'''RUN_TIRSHAKER: The Tirshaker mode {cfg.tirshaker.mode} is not yet fully functional. Please use a different mode
''')
    
    

    
    if cfg.general.cube_name == None:
        #We assume the cube path is in inset
        Tirific_Template['INSET'] = f"../{Tirific_Template['INSET']}"
    else:
        Tirific_Template['INSET'] = f"{cfg.general.cube_name}"

    if cfg.general.multiprocessing:
        processes = 1
        tmpcpu = cfg.general.ncpu
        tmpcpu -= int(float(Tirific_Template['NCORES']))
        while  tmpcpu  >= int(float(Tirific_Template['NCORES'])):
            processes += 1
            tmpcpu  -= int(float(Tirific_Template['NCORES']))
        if processes == 1:
           Tirific_Template['NCORES'] = f'{cfg.general.ncpu}' 

    else: 
        processes = 1
        if cfg.general.ncpu != -1:
            if cfg.general.ncpu < 10:
                Tirific_Template['NCORES'] = cfg.general.ncpu
            else:
                Tirific_Template['NCORES'] = 10
  
    Tirific_Template['OUTSET'] = ''
    Tirific_Template['PROGRESSLOG'] = ''
    Tirific_Template['TEXTLOG'] = ''
    Tirific_Template['TIRSMO'] = ''
    Tirific_Template['COOLGAL'] = ''
    Tirific_Template['TILT'] = ''
    Tirific_Template['BIGTILT'] = ''
    if processes == 1:
        if cfg.tirshaker.inimode != -1:
            Tirific_Template['INIMODE'] = cfg.tirshaker.inimode
        else:
            if Tirific_Template['INIMODE'] == '':
                Tirific_Template['INIMODE'] = 0 
    else:
        Tirific_Template['INIMODE'] = 0

    Tirific_Template['LOGNAME'] = 'Error_Shaker.log'
    Tirific_Template['TIRDEF'] = 'Error_Shaker_Out.def'
    out = [f'Parameter = {x} with block = {fit_groups[x]["BLOCK"]} for the rings {fit_groups[x]["RINGS"]} and disks {fit_groups[x]["DISKS"]} varied by {fit_groups[x]["VARIATION"][0]}. \n' for x in fit_groups ]
    if cfg.general.verbose:
        log_statement += print_log(f'''RUN_TIRSHAKER: We are shaking with the following parameters:
{''.join(out)}
''',log)
    fit_groups['TO_COLLECT'] = []
    fit_groups['COLLECTED'] = {}
    for group in  fit_groups:
        if group not in ['COLLECTED','TO_COLLECT']:
            for disk in fit_groups[group]['DISKS']:
                para = group.split('_')[0]
                if disk != 1:
                    para = f'{para}_{disk}'
                if para not in  fit_groups['TO_COLLECT']:
                    fit_groups['TO_COLLECT'].append(para)
                fit_groups['COLLECTED'][para] = []
    
    return log_statement,Tirific_Template,fit_groups,processes
