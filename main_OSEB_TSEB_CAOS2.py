#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run TSEB/OSEB over a timeseries of point measurements
for the data set of Petit-Nobressart from summer 2015 collected during 
the CAOS field campaign. This repository is part of the HESS submission
"Understanding biases in evapotranspiration schemes using a pattern-based 
evaluation of the diurnal cycle: a case study in Luxembourg" 2018.

To run the models two files have to be provided:
- Configuration file including location data, emissivity, spectral data,
     information on G computation method etc. The configuration file is a
     dictionary containing the required data. The file should be named
     setup_config.py.
     
- Timeseries of meteorological and vegetation properties.
     For the OSEB/TSEB model the following variables are required:
     u (wind speed), Ta (air temperature), rH (relative humidity),
     p (atmospheric pressure), Sdn (shortwave downwelling radiation),
     Ldn (longwave downwelling radiatio), es (saturation vapor 
     pressure deficit), G (soil heat flux), Trad (radiometric surface 
     temperature).
     In the case that vegetation properties are time-dependent they also
     have to be specified in this file, otherwiese they can be set as
     scalar values in the configuration file.

The data to run the timeseries for the CAOS 2 field campaign is not included
in this repository, but has to be requested from the corresponding author 
of the related publication.

"""

import argparse
import sys
from pathlib import Path

HOME = Path(__file__).parent
sys.path.append(str(HOME / 'pyTSEB_fork'))
from src.pyTSEB import PyTSEB
    
def get_args():
    ''' Parse input file name, model name and configuration file name from
    command line input arguments.
    
    If no arguments are parsed then the values in the configuration file are
    used.
    
    '''
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input_file', 
                        help=('File with meteorological and vegetation ' + \
                        'properties'), required=False)
    parser.add_argument('-m', '--model', help='Choose either TSEB or OSEB',
                        required=False)
    
    args = parser.parse_args()
    input_file = args.input_file
    model_name = args.model
    
    if input_file is not None:
        if not (HOME / 'model_input' / input_file).is_file():
            raise FileNotFoundError(f'Could not find file at {input_file}.')
        input_file = Path(input_file)
        
    if model_name is not None:
        if (model_name != 'TSEB') and (model_name != 'OSEB'):
            raise ValueError('For the model choose either TSEB or OSEB.')
    
    return input_file, model_name
    
    
def runTSEB_from_config_file(config_file):
    ''' Runs the TSEB/OSEB model'''
    # Create a class instance from PyTSEB
    setup = PyTSEB()
    # Write configuration data into PyTSEB instance
    setup.GetDataTSEB(config_file, isImage=False)
    # Run model in timeseries mode
    setup.RunTSEBPointSeriesArray()
    return


def main(input_file, model_name):     
    # Loading config file
    sys.path.append(str(HOME))
    from model_input.setup_config import CONFIG_FILE 
    
    if model_name is None:
        model_name = CONFIG_FILE['ModelMode'][:4]
    if input_file is None:
        input_file = CONFIG_FILE['PointTimeseriesInput']['InputFile']
        
    output_name = f"model_output_timeseries_{model_name}.txt"   
    
    # Write correct input name into the configuration file
    CONFIG_FILE['PointTimeseriesInput']['InputFile'] = (HOME /
               'model_input' / input_file)
    # Write correct output name into the configuration file
    CONFIG_FILE['PointTimeseriesInput']['OutputFile'] = (HOME /
               'model_output' / output_name)
    
    CONFIG_FILE['ModelMode'] = model_name
    if model_name == 'TSEB':
        CONFIG_FILE['ModelMode'] = 'TSEB_PT'
        
    # Create output directory if it does not exist
    if not (HOME / 'model_output').exists():
        (HOME / 'model_output').mkdir()
    
    print(model_name)
    print(input_file)
    # Run the model
    runTSEB_from_config_file(CONFIG_FILE)
    
if __name__ == '__main__':
    input_file, model_name = get_args()
    main(input_file, model_name)



