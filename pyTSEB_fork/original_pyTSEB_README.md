# PyTSEB

## Synopsis

This project contains *Python* code for *Two Source Energy Balance* models (Priestley-Taylor **TSEB-PT**, 
Dual Time Difference **DTD** and TSEB with component soil and canopy temperatures **TSEB-2T**) 
for estimating sensible and latent heat flux (evapotranspiration) based on measurements of radiometric surface temperature. 

The project consists of: 

1. lower-level modules with the basic functions needed in any resistance energy balance model 

2. higher-level scripts for easily running TSEB with tabulated data and/or satellite/airborne imagery.

## Installation

Download the project to your local system, enter the download directory and then type

`python setup.py install` 

if you want to install pyTSEB and its low-level modules in your Python distribution. 

The following Python libraries will be required for running TSEB over an image:

- Numpy
- GDAL

## Code Example
### High-level example

The easiest way to get a feeling of TSEB and its configuration is through the provided ipython/jupyter notebooks. 
In a terminal shell, navigate to your working folder and type

- `jupyter notebook ProcessPointTimeSeries.ipynb` 
>for configuring and running TSEB over a time series of tabulated data

- `jupyter notebook ProcessLocalImage.ipynb` 
>for configuring and running TSEB over an image/scene using local meteorological data

In addition, you can also run TSEB with the scripts *MAIN_TSEB_LocalImage.py* and *MAIN_TSEB_PointTimeSeries.py*, 
which will read an input configuration file (defaults are *Config_LocalImage.txt* and *Config_PointTimeSeries.txt* respectively). 
You can edit these configuration files or make a copy to fit your data and site characteristics and either run any of 
these two scripts in a Python GUI or in a terminal shell:

- `python MAIN_TSEB_LocalImage.py <configuration file>`
> where \<configuration file> points to a customized configuration file... leave it blank if you want to use the default 
file *Config_LocalImage.txt*

- `python MAIN_TSEB_PointTimeSeries.py <configuration file>`
> where \<configuration file> points to a customized configuration file... leave it blank if you want to use the default 
file *Config_PointTimeSeries.txt*

### Low-level example
You can run any TSEB model or any related process in python by importing the module *TSEB*. 
It will also import the ancillary modules (*resitances.py* as `res`, *netRadiation* as `rad`,
*MOsimilarity.py* as `MO`, *ClumpingIndex.py* as `CI` and *meteoUtils.py* as `met`)

```python
import TSEB 
output=TSEB.TSEB_PT(Tr_K, vza, Ta_K, u, ea, p, Sdn_dir, Sdn_dif, fvis, fnir, sza, Lsky, LAI, hc, emisVeg, emisGrd, spectraVeg, spectraGrd, z_0M, d_0, zu, zt)
```

You can type
`help(TSEB.TSEB_PT)`
to understand better the inputs needed and the outputs returned

The direct and difuse shortwave radiation (`Sdn_dir`, `Sdn_dif`, `fvis`, `fnir`) and the downwelling longwave radiation (`Lsky`) can be estimated by

```python
emisAtm = TSEB.rad.CalcEmiss_atm(ea,Ta_K_1) # Estimate atmospheric emissivity from vapour pressure (mb) and air Temperature (K)
Lsky = emisAtm * TSEB.met.CalcStephanBoltzmann(Ta_K_1) # in W m-2
difvis,difnir, fvis,fnir=TSEB.rad.CalcDifuseRatio(Sdn,sza,press=p, Wv=1) # fraction of difuse and PAR/NIR radiation from shortwave irradiance (W m-2, solar zenith angle, atmospheric pressure and precipitable water vapour )
Skyl=difvis*fvis+difnir*fnir # broadband difuse fraction
Sdn_dir=Sdn*(1.0-Skyl)
Sdn_dif=Sdn*Skyl
```
   
## Basic Contents
### High-level modules
- *.src/pyTSEB.py*, class object for TSEB scripting

- *ProcessPointTimeSeries.ipynb* and *ProcessLocalImage.ipynb* notebooks for using TSEB and configuring 
TSEB through a Graphical User Interface, GUI

- *MAIN_TSEB_LocalImage.py* and *MAIN_TSEB_PointTimeSeries.py*, high level scripts for running TSEB 
through a configuration file (*Config_LocalImage.txt* or *Config_PointTimeSeries.txt*)

- *TSEB_ProcessPointTimeSeries.ipynb* notebook with more details about the high-level TSEB programming

- *pyTSEB_in_Detail.ipynb* notebook with a closer look at the low-level and high-level TSEB code

- *TSEB_and_Resistances.ipynb* notebook for undertanding the TSEB model and the estimation of resistances

### Low-level modules
The low-level modules in this project are aimed at providing customisation and more flexibility in running TSEB. 
The following modules are included

- *.src/TSEB.py*
> core functions for running different TSEB models (`TSEB_PT (*args,**kwargs)`, `TSEB_2T(*args,**kwargs)`, 
`DTD (*args,**kwargs)`), or a One Source Energy Balance model (`OSEB(*args,**kwargs)`). 

- *.src/netRadiation.py*
> functions for estimating net radiation and its partitioning between soil and canopy

- *.src/resistances.py*
> functions for estimating the different resistances for momemtum and heat transport and surface roughness

- *.src/MOsimilarity.py*
> functions for computing adiabatic corrections for heat and momentum transport, 
Monin-Obukhov length, friction velocity and wind profiles

- *.src/ClumpingIndex.py*
> functions for estimating the canopy clumping index and get effective values of Leaf Area Index

- *.src/meteoUtils.py*
> functions for estimating meteorolgical-related variables such as density of air, 
heat capacity of air or latent heat of vaporization.

## API Reference
http://pytseb.readthedocs.org/en/latest/index.html

## Main Scientific References
- Norman,  J.  M.,  Kustas,  W.  P.,  Prueger,  J.  H.,  and  Diak,  G.  R.: Surface  flux  estimation  using  radiometric  temperature:  a  dual-temperature-difference method to minimize measurement errors, Water  Resour.  Res.,  36,  2263,  doi: 10.1029/2000WR900033, 2000
- Norman,  J.,  Kustas,  W.,  and  Humes,  K.:  A  two-source  approach for estimating soil and vegetation fluxes from observations of directional radiometric surface temperature, Agr. Forest Meteorol., 77, 263–293, doi: 10.1016/0168-1923(95)02265-Y, 1995
- Kustas, W. P. and Norman, J. M.: A two-source approach for estimating turbulent fluxes using multiple angle thermal infrared observations, Water Resour. Res., 33, 1495–1508, 199
- Kustas,  W.  P.  and  Norman,  J.  M.:  Evaluation  of  soil  and  vegetation heat flux prediction using a simple two-source model with radiometric  temperatures  for  partial  canopy  cover,  Agr.  Forest Meteorol., 94, 13–29, 199

## Tests
The folder *./Input* contains examples for running TSEB in a tabulated time series (*ExampleTableInput.txt*) 
and in an image (*ExampleImage_\< variable >.tif*). Just run the high-level scripts with the configuration files 
provided by default and compare the resulting outputs with the files stored in *./Output/*

## Contributors
- **Hector Nieto** <hnieto@ias.csic.es> <hector.nieto.solana@gmail.com> main developer
- **William P. Kustas** TSEB modeling, tester 
- **Radoslaw Guzinski** DTD code developer, tester
- **Ana Andreu** tester

## License
pyTSEB: a Python Two Source Energy Balance Model

Copyright 2016 Hector Nieto and contributors.
    
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
