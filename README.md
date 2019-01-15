# ***PyTSEB model for the CAOS 2 field campaign***

## Overview

This repository provides the code to run the One- or Two-Source Energy Balance (OSEB/TSEB) model used in the 2018 HESS submission "Understanding biases in evapotranspiration schemes using a pattern-based evaluation of the diurnal cycle: a case study in Luxembourg" by Maik Renner (MPI-BGC), Claire Brenner (BOKU) , Kaniska Mallick (LIST), Hans-Dieter Wizemann (UHOH), Luigi Conte (MPI-BGC), Ivonne Trebs (LIST), Jianhui Wei (UHOH),  Volker Wulfmeyer (UHOH), Karsten Schulz (BOKU), Axel Kleidon (MPI-BGC).

## pyTSEB

The main body of the code is a fork of the pyTSEB repository of hectornieto from this
[commit](https://github.com/hectornieto/pyTSEB/tree/ac3fe785ead3a9c4f09d773e00a9334f75f21fd0). The copyright for the code is with Hector Nieto 2016. Further information on the TSEB/OSEB model as implemented in pyTSEB is provided [here](http://pytseb.readthedocs.org/en/latest/index.html). The original README is included at pyTSEB_fork/original_pyTSEB_README.md.

## Model run options

To run the main file main_OSEB_TSEB_CAOS2.py two inputs have to be provided:
* ***Configuration file*** including location data, emissivity, spectral data,
     information on G computation method etc. The configuration file is a
     dictionary containing the required data. The file should be named
     setup_config.py and saved in the folder model_input.
* ***Timeseries of meteorological and vegetation properties*** with the following variables: <br>
u (wind speed), Ta (air temperature), rH (relative humidity),
p (atmospheric pressure), Sdn (shortwave downwelling radiation),
Ldn (longwave downwelling radiatio), es (saturation vapor
pressure deficit), G (soil heat flux), Trad (radiometric surface
temperature). <br>
rH and es are used to calculate ea (atmospheric vapor pressure),
which is used to calculate cp (heat capacity of dry air) and
rho (air density).

The data to run the timeseries for the CAOS 2 field campaign is not included
in this repository, but has to be requested from the corresponding author
of the related HESS publication. Templates with the required variables and their labels for both input files are given in the folder model_input.


The main file can be run from the command line using the following arguments_
-f meteorological input file name (which should be saved in the folder model_input)
-m model model (either OSEB or TSEB)
In the case that no arguments are provided the values from the configuration file are used.

[![DOI](https://zenodo.org/badge/134852988.svg)](https://zenodo.org/badge/latestdoi/134852988)


