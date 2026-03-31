from dataclasses import dataclass, field
import psutil
from omegaconf import MISSING, OmegaConf
from typing import List, Optional
import os
import sys
from TRM_errors.config.variables_config import Min_Errors,Variations
#The default total database currently make 229 galaxies
@dataclass
class Tirshaker:
    #tirshaker is the default
    enable: bool = True
    #The input def file for which to calculate the errors
    deffile_in: str = 'Finalmodel.def'
    deffile_out: str = 'Shaken_Errors.def'
    directory: str = 'Error_Shaker'
    #Do we want a log
    log: bool = False 
    mode: str = 'fitted' #Fitted the settings and grouping will be read from the fits file and def file if manual they have to be provided
    inimode: int=-1
    iterations: int=20
    individual_loops: int = -1  #Set this to -1 for final release
    tirific: str = 'tirific'

@dataclass
class General:
    cube_name: Optional[str] = None
    verbose: bool = True
    try:
        ncpu: int = len(psutil.Process().cpu_affinity())-1
    except AttributeError:
        ncpu: int = psutil.cpu_count()-1
   
    directory: str = os.getcwd()
    multiprocessing: bool = True
    calc_mode: str = 'mad'
    clean: bool = True
    #font_file: str = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"
@dataclass
class defaults:
    print_examples: bool = False
    configuration_file: Optional[str] = None
    general: General =  field(default_factory = General)
    tirshaker: Tirshaker=  field(default_factory = Tirshaker)
    min_errors: Min_Errors =  field(default_factory = Min_Errors)
    variations: Variations =  field(default_factory = Variations) 


def load_config(argv):
    cfg = OmegaConf.structured(defaults)
    if cfg.general.ncpu == psutil.cpu_count():
        cfg.general.ncpu -= 1
    inputconf = OmegaConf.from_cli(argv)
    cfg_input = OmegaConf.merge(cfg,inputconf)
    
    if cfg_input.print_examples:
        with open('TRM_errors_default.yml','w') as default_write:
            default_write.write(OmegaConf.to_yaml(cfg))
        print(f'''We have printed the file TRM_errors_default.yml in {os.getcwd()}.
Exiting moments.''')
        sys.exit()

    if cfg_input.configuration_file:
        succes = False
        while not succes:
            try:
                yaml_config = OmegaConf.load(cfg_input.configuration_file)
        #merge yml file with defaults
                cfg = OmegaConf.merge(cfg,yaml_config)
                succes = True
            except FileNotFoundError:
                cfg_input.configuration_file = input(f'''
You have provided a config file ({cfg_input.configuration_file}) but it can't be found.
If you want to provide a config file please give the correct name.
Else press CTRL-C to abort.
configuration_file = ''')
    cfg = OmegaConf.merge(cfg,inputconf) 


    if not cfg.general.cube_name is None:
        cfg = check_cube_name(cfg)
       
    return cfg

def check_cube_name(cfg):
    directory = os.path.dirname(cfg.general.cube_name)
    if directory == '':
        cfg.general.cube_name = os.path.join(cfg.general.directory,
            cfg.general.cube_name)
    return cfg