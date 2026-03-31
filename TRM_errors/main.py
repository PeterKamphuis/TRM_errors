# -*- coding: future_fstrings -*-

# This is the stand alone version of the pyFAT moments to create moment maps

#from optparse import OptionParser

from TRM_errors.config.config import load_config
from TRM_errors.common.common import load_tirific,check_cpu
import numpy as np
import sys
import os
import traceback
import warnings
import TRM_errors
import psutil
from multiprocessing import get_context,Manager

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

def main_prog():
    main(sys.argv[1:])



def main(argv):
    if '-v' in argv or '--version' in argv:
        print(f"This is version {TRM_errors.__version__} of the program.")
        sys.exit()

    if '-h' in argv or '--help' in argv:
        print('''
Use TRM_Errors in this way:

All config parameters can be set directly from the command line by setting the correct parameters, e.g:
              
create_TRM_errors def_file=cube.fits error_generator=tirshaker 

              ''')
        sys.exit()

    cfg = load_config(argv)

    cfg = check_cpu(cfg)
    
    # for some dumb reason pools have to be called from main
    if cfg.tirshaker.enable:
        from TRM_errors.tirshaker.tirshaker import prepare_template, set_individual_iteration,\
                    run_tirific,run_individual_iteration,tirshaker_cleanup,finish_current_run
        #First prepare the main template
        log_statement,Tirific_Template,fit_groups, no_processes = prepare_template(cfg)

        if no_processes == 1:
            current_run = 'not_set'
            # start loop for all iterations
            for i in range(cfg.tirshaker.iterations):

                out = set_individual_iteration(Tirific_Template, i,fit_groups,\
                            f'{cfg.general.directory}/{cfg.tirshaker.directory}', cfg.tirshaker.tirific,\
                            name_in=f'Error_Shaker_In.def',verbose=cfg.general.verbose,clean= cfg.general.clean)
                log_statement += out['log']
               
                current_run = run_tirific(current_run,deffile=out['deffile'],work_dir = out['directory']\
                                    ,tirific_call= cfg.tirshaker.tirific, \
                                    max_ini_time= int(300*(int(Tirific_Template['INIMODE'])+1)))
                if not cfg.general.clean:
                    source_in = f"{out['directory']}/{out['deffile']}"
                    target_in = f"{out['directory']}/Error_Shaker_In_{i}.def"
                    source_out = f"{out['directory']}/{out['tmp_name_out']}"
                    target_out = f"{out['directory']}/Error_Shaker_Out_{i}.def"
                    os.system(f'''cp {source_in} {target_in}''')
                    os.system(f'''cp {source_out} {target_out}''')

            
            # Read the values of the pararameter groups
                for parameter in fit_groups['TO_COLLECT']:
                    fit_groups['COLLECTED'][parameter].append(load_tirific(f"{out['directory']}/{out['tmp_name_out']}",\
                        Variables = [parameter],array=True))
            log_statement += finish_current_run(current_run)

        else:
             #In case of multiprocessing we want to make sure to start with
        #The big galaxies
        #Setup an array of configs with locks
            list_iterations = []
            for i in range(cfg.tirshaker.iterations):
                list_iterations.append([set_individual_iteration(Tirific_Template, i,fit_groups, \
                                    f'{cfg.general.directory}/{cfg.tirshaker.directory}', \
                                    cfg.tirshaker.tirific,name_in=f'Error_Shaker_In_{i}.def',\
                                    name_out=f'Error_Shaker_Out_{i}.def',
                                    verbose=cfg.general.verbose,clean= cfg.general.clean) \
                                    ])
            # simply running tiriic would mess up with continuation/restart id
            with Manager() as loop_manager:
                with get_context("spawn").Pool(processes=no_processes) as pool:
                    if cfg.general.verbose:
                        print(f'Starting iterations with {no_processes} processes')
                    iteration_output = pool.starmap(run_individual_iteration, list_iterations)
            for iter in iteration_output:
                log_statement += iter['log']
                for parameter in fit_groups['TO_COLLECT']:
                    fit_groups['COLLECTED'][parameter].append(iter[parameter])
        # and then finish up
        log_statement += tirshaker_cleanup(fit_groups, cfg)
        







if __name__ =="__main__":
    main()
