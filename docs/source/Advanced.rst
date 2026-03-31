.. TRM_errors documentation master file, created by
   sphinx-quickstart on Tue Oct 31 17:10:19 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to TRM_errors's documentation!
======================================


TRM_errors errors is a package to create actual errors for Tilted ring models.

The options to create models are 

Tirshaker
----------------
In order to create errors run create_TRM_errors.

It takes the following command line or yaml input

as singular input commands

**print_examples**
    *bool = False*

    print an example yaml file

**configuration_file**
    
    *Optional[str] = None*

    Input configuration yaml file

Tirshaker input 
------------------------------------

*(tirshaker.command=)*

**enable**

    *bool = True*

    Run true or False  (Tirshaker is the default error code)


**deffile_in**

    *str = 'Finalmodel.def'*

    The input def file for which to calculate the errors


**deffile_out**
    
    *str = 'Shaken_Errors.def'*
    
**directory**

    *str = 'Error_Shaker'*

**log**

    *bool = False* 

**mode**

    *str = 'fitted'*
     
     'fitted'= the settings and grouping will be read from the def file.
      !!!!!!!!!!!  Manual is not functioning yet as the fitting parameters are not adapted!!!!!!!!!
     'manual' = The variations and groups are read from the yaml file. 
       !!!!!!!!!!!  Manual is not functioning yet as the fitting parameters are not adapted!!!!!!!!!  

**inimode**

    *int=2*

    initiazation mode of runs in manual mode


**iterations**

    *int=20*

    Amount of shaker iterations

**individual_loops**

    *int = 3*
    
    Amount of individual loops in each iteration

**tirific**

    *str = 'tirific'*

    command to run tirific


General input 
---------------------------------
*(general.command=)*

    **cube_name**
    
     *Optional[str] = None*

    The name of the cube to be used in the input def file. If not provided it taken from the input def file. If a single file name it is assumed to be in the main directory.

    **ncpu**
    
    *int = cpu_count()-1*

    Number of cpu to use for multiprocessing, default is all but one. If set to 1 then multiprocessing is turned off.
    
    **directory**
    
    *str = os.getcwd()*

    Main directory where the input def file is located and where the output will be stored. If not provided it is assumed to be the current working directory.
    
    **verbose**
    
    *bool = True*

    print info or not

    **calc_mode**

    *str =  'mad'

    analysis of the output, mad provides a filtered std over the iterations, 'fat' includes the difference beteen the final output and the median output in the itreations into the calculation. 

    **multiprocessing**
    
    *bool = True*
    
    Use multiple cpus

    **clean**

    *bool = True*

    If set to false then all input and output defiles will be maintain in the output directory


Variations input
-------------------------------------
*(variation.command=)*

    **VARY**
    
    *str = '!VROT VROT_2'*    
     
    Set the parameters manually following tirific VARY syntax, if ! causes problems use i 
    
   
    The parameters can set as a multiple of the resolution or as  absolute value in tirific units.
    They are a list of 4 variables where the first indicates the variation,  the second whether this is a multiple of the resolution (res) or an absoluter value (unit). Then the unit of the parameter.
    If the unit is angle there is no associated resolution. For km/s or m/s the resolution is the channel width and for  degree, armin or arcsec we use the beam.
    for 'jy/arcsec^2' the noise is considered the resolution. The final parameter should for now always be 'a'.

   
    **PA**
    
    *List = field(default_factory=lambda: [10, 'unit','angle','a'])*
    
    **INCL**
    
    *List = field(default_factory=lambda: [10, 'unit','angle', 'a'])*
    
    **VROT**
    
    *List = field(default_factory=lambda: [5, 'res','km/s','a'])*
    
    **VRAD**
    
    *List = field(default_factory=lambda: [2.5, 'res','km/s','a'])*    
    
    **VSYS**
    
    *List = field(default_factory=lambda: [0.1, 'res','km/s', 'a'])*
    
    **XPOS**
    
    *List = field(default_factory=lambda: [0.3, 'res','degree','a'])*
    
    **YPOS**
    
    *List = field(default_factory=lambda: [0.3, 'res','degree','a'])*    
    
    **SBR**
    
    *List = field(default_factory=lambda: [1e-4, 'unit','jy/arcsec^2', 'a'])*
    
    **Z0**
    
    *List = field(default_factory=lambda: [1, 'res','arcsec','a'])*
    
    **SDIS**
    
    *List = field(default_factory=lambda: [2, 'res','km/s','a'])*    

minimum errors
-------------------------------------
*(min_errors.parameter =)*

    **PA** 
    
    *float = 0.*
    
    **INCL**

    *float = 0.*
    
    **VROT**

    *float = 0.*
    
    **VRAD**
    
    *float = 0.*
    
    **VSYS**
    
    *float = 0.*
    
    **XPOS**
    
    *float = 0.*
    
    **YPOS**
    
    *float = 0.*
    
    **SBR** 
    
    *float = 0.*
    
    **Z0** 
    
    *float = 0.*
    
    **SDIS** 
    
    *float = 0.*   

.. toctree::
   :maxdepth: 2
   :caption: Contents:


