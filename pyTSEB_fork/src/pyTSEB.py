# This file is part PyTSEB, consisting of of high level pyTSEB scripting
# Copyright 2016 Hector Nieto and contributors listed in the README.md file.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Created on Thu Jan  7 16:37:45 2016
@author: Hector Nieto (hnieto@ias.csic.es)

Modified on Jan  27 2016
@author: Hector Nieto (hnieto@ias.csic.es)

DESCRIPTION
===========
This package contains the class object for configuring and running TSEB for both
an image with constant meteorology forcing and a time-series of tabulated data.

EXAMPLES
========
The easiest way to get a feeling of TSEB and its configuration is throuh the ipython/jupyter notebooks. 

Jupyter notebook pyTSEB GUI
---------------------------
To configure TSEB for processing a time series of tabulated data, type in a ipython terminal or a jupyter notebook.

.. code-block:: ipython

    import pyTSEB # Import pyTSEB module
    tseb=pyTSEB.PyTSEB() # Create a PyTSEB instance
    tseb.PointTimeSeriesWidget() # Launches the GUI

then to run pyTSEB.

.. code-block:: ipython

    tseb.GetDataTSEBWidgets(isImage=False)# Get the data from the widgets
    tseb.RunTSEBLocalImage()# Run TSEB

Similarly, to configure and run TSEB for an image.

.. code-block:: ipython

    import pyTSEB # Import pyTSEB module
    tseb=pyTSEB.PyTSEB() # Create a PyTSEB instance
    tseb.PointLocalImageWidget() # Launches the GUI
    tseb.GetDataTSEBWidgets(isImage=True)# Get the data from the widgets
    tseb.RunTSEBLocalImage()# Run TSEB

Parsing directly a configuration file
-------------------------------------
You can also parse direcly into TSEB a configuration file previouly created.

>>> tseb=PyTSEB()
>>> configData=tseb.parseInputConfig(configFile,isImage=True) #Read the data from the configuration file into a python dictionary
>>> tseb.GetDataTSEB(configData,isImage=True) # Parse the data from the dictionary to TSEB
>>> tseb.RunTSEBLocalImage() # Run TSEB

see the guidelines for input and configuration file preparation in :doc:`README_Notebooks`.

"""

class PyTSEB():
    
    def __init__(self):
        '''Initialize input variables  with default  values'''
        self.InputFile='./Input/ExampleInput.txt'
        self.OutputTxtFile='./Output/test.txt'
        self.OutputImageFile='./Output/test.tif'
        
        # MOdel to run
        self.model='TSEB_PT'
        # Site Description
        self.lat=36.95
        self.lon=2.33
        self.alt=200
        self.stdlon=15
        self.zu=2.5
        self.zt=3.5
        # Spectral Properties
        self.rhovis=0.07
        self.rhonir=0.32
        self.tauvis=0.08
        self.taunir=0.33
        self.rsoilvis=0.15
        self.rsoilnir=0.25
        self.emisVeg=0.98
        self.emisGrd=0.95
        # Surface Properties
        self.max_PT=1.26
        self.x_LAD=1.0
        self.leaf_width=0.1
        self.z0soil=0.01
        self.LANDCOVER=11
        # Soil Heat Flux calculation
        self.CalcG=1
        self.Gconstant=0
        self.Gratio=0.35
        self.GAmp=0.35
        self.Gphase=3
        self.Gshape=24
        # Default Vegetation variables
        self.f_c=1.0
        self.f_g=1.0
        self.wc=1.0
        # Output variables saved in images
        self.fields=('H1','LE1','R_n1','G1')
        # Ancillary output variables
        self.anc_fields=('H_C1','LE_C1','LE_partition','T_C1', 'T_S1','R_ns1','R_nl1', 'u_friction', 'L', 'R_A1', 'R_X1', 'R_S1', 'flag')
        # File Configuration variables
        self.input_commom_vars=('TSEB_MODEL','lat','lon','altitude','stdlon',
                   'z_t','z_u','emisVeg','emisGrd','rhovis','tauvis',
                   'rhonir','taunir','rsoilvis','rsoilnir','Max_alpha_PT',
                   'x_LAD','z0_soil','LANDCOVER','leaf_width',
                   'CalcG','G_constant','G_ratio','GAmp','Gphase','Gshape','OutputFile')
        self.input_image_vars=['Input_LST','Input_VZA','USE_MASK','Input_LAI','Input_Fc',
                         'Input_Fg','Input_Hc','Input_Wc','OutputFile','Time','DOY','Ta_1',
                         'Ta_0','u','ea', 'Sdn', 'Ldn','p']
        self.input_point_vars=['InputFile']
        self.resolution_dependence=False
           
    def GetDataTSEB(self, config_file, isImage=True):
        '''
        Parses the parameters in a configuration file directly to 
        TSEB variables for running TSEB
        '''
        self.TSEB_MODEL=config_file['ModelMode']
        
        # Site properties
        SiteProps = config_file['SiteProps']
        self.lat,self.lon,self.alt,self.stdlon,self.zu,self.zt = (float(SiteProps['lat']),
                float(SiteProps['lon']),float(SiteProps['altitude']),float(SiteProps['stdlon']),
                float(SiteProps['zu']),float(SiteProps['zt']))
        
        # Emissivities
        SiteVegProps = config_file['SiteVegProps']
        self.emisVeg,self.emisSoil = (float(SiteVegProps['emisVeg']), 
                float(SiteVegProps['emisSoil']))
        # Canopy spectral properties   
        self.spectraVeg = {'rho_leaf_vis':float(SiteVegProps['rhovis']), 'tau_leaf_vis':float(SiteVegProps['tauvis']),
                    'rho_leaf_nir':float(SiteVegProps['rhonir']), 'tau_leaf_nir':float(SiteVegProps['taunir'])}
        # Soil spectral properties
        self.spectraGrd = {'rsoilvis':float(SiteVegProps['rsoilvis']), 'rsoilnir':float(SiteVegProps['rsoilnir'])}
        
        # Ancillary data
        AncProps = config_file['AncProps']
        self.Max_alpha_PT, self.x_LAD, self.leaf_width,self.z0_soil,self.LANDCOVER = (float(AncProps['maxAlphaPT']),
                        float(SiteVegProps['xLAD']),float(SiteVegProps['leafWidth']),float(AncProps['z0soil']),int(SiteProps['landCover']))
        
        # Soil heat flux paramerization
        SHFinfo = config_file['SHFinfo']
        if int(SHFinfo['CalcG'])==0:
            self.CalcG=[0,float(SHFinfo['G_constant'])]
        elif int(SHFinfo['CalcG'])==1:
            self.CalcG=[1,float(SHFinfo['G_ratio'])]
        elif int(SHFinfo['CalcG'])==2:
            # The twelve has to be overwritten
            self.CalcG=[2,[12.0,float(SHFinfo['GAmp']),float(SHFinfo['Gphase']),float(SHFinfo['Gshape'])]]
        
        if isImage:
            Imagery = config_file['Imagery']
            self.OutputFile=Imagery['outputFile']
            # Get all the input parameters
            self.input_LST=str(Imagery['inputLST']).strip('"')
            self.input_VZA=str(AncProps['VZA']).strip('"')
            self.input_LAI=str(SiteVegProps['LAI']).strip('"')
            self.input_Hc=str(SiteVegProps['Hc']).strip('"')
            self.input_Fc=str(SiteVegProps['Fc']).strip('"')
            self.input_Fg=str(SiteVegProps['Fg']).strip('"')
            self.input_Wc=str(SiteVegProps['Wc']).strip('"')
            self.input_mask=str(AncProps['useMask']).strip('"')
            
            MeteoData = config_file['Meteo']
            self.DOY,self.Time,self.Ta_1,self.Sdn,self.u,self.ea,self.Ldn,self.p=(float(MeteoData['DOY']),float(MeteoData['Time']),
                float(MeteoData['Ta_1']),float(MeteoData['Sdn']),float(MeteoData['u']),float(MeteoData['ea']),
                str(MeteoData['Ldn']).strip('"'),str(MeteoData['p']).strip('"'))
            
            # Get correct time for time shift according to Santanello and Friedl
            if int(SHFinfo['CalcG'])==2:
                # The twelve has to be overwritten
                self.CalcG[1][0] = self.Time
            
            if self.TSEB_MODEL=='DTD':
                self.Ta_0=float(MeteoData['Ta_0'])
            
            modelMode = self.TSEB_MODEL
            if modelMode == 'OSEB' or modelMode == 'DTD_OSEB':
                self.input_LAI = str(0.0)
                if modelMode == 'DTD_OSEB':
                    modelMode = 'DTD'
                    self.TSEB_MODEL = 'DTD'
                else:
                    modelMode = 'TSEB_PT'
                    self.TSEB_MODEL = 'TSEB_PT'  
                if int(SHFinfo['CalcG'])==1:
                    self.CalcG=[1,float(0.15)]
                    
            self.resolution_dependence = config_file['resolution_dep']
            
        else:
            PointTimeseriesInput = config_file['PointTimeseriesInput']
            self.InputFile=str(PointTimeseriesInput['InputFile']).strip('"')
            self.OutputFile=str(PointTimeseriesInput['OutputFile']).strip('"')
            self.f_c = float(SiteVegProps['Fc'])
            self.f_g = float(SiteVegProps['Fg'])
            self.wc = float(SiteVegProps['Wc'])   
            self.kB_method = PointTimeseriesInput['kB']
    
    def CheckDataPointSeriers(self):
        '''Checks that all the data required for TSEB is contained in an input ascci table'''
        success=False
        # Mandatory Input Fields
        MandatoryFields_TSEB_PT=('Year','DOY','Time','Trad','VZA','Ta','u','ea','Sdn','LAI','hc')
        MandatoryFields_DTD=('Year','DOY','Time','Trad_0','Trad','VZA','Ta_0','Ta','u','ea','Sdn','LAI','hc')                        
        MandatoryFields_TSEB_2T=('Year','DOY','Time','Tc','Ts','Ta','u','ea','Sdn','LAI','hc')  
        # Check that all mandatory input variables exist
        if self.TSEB_MODEL=='TSEB_PT':
            for field in MandatoryFields_TSEB_PT:
                if field not in self.inputNames:
                    print('ERROR: ' +field +' not found in file '+ self.InputFile)
                    return success
        elif self.TSEB_MODEL=='OSEB':
            for field in MandatoryFields_TSEB_PT:
                if field not in self.inputNames:
                    print('ERROR: ' +field +' not found in file '+ self.InputFile)
                    return success
        elif self.TSEB_MODEL=='DTD':
            for field in MandatoryFields_DTD:
                if field not in self.inputNames:
                    print('ERROR: ' +field +' not found in file '+ self.InputFile)
                    return success
        elif self.TSEB_MODEL=='TSEB_2T':
            for field in MandatoryFields_TSEB_2T:
                if field not in self.inputNames:
                    print('ERROR: ' +field +' not found in file '+ self.InputFile)
                    return success
        else:
            print('Not valid TSEB model, check your configuration file')
            return success
        return True

    def RunTSEBPointSeriesArray(self):
        ''' Runs TSEB for all the pixel in an image'''
        import src.TSEB as TSEB
        from  os.path import dirname, exists
        from os import mkdir
        import numpy as np  
        import csv
        
        def addData(dataDict, fieldName, fieldValue):
            if fieldName not in dataDict.keys():
                dataDict[fieldName] = np.array([])
            dataDict[fieldName] = np.append(dataDict[fieldName], fieldValue)


        #======================================
        # Process input file

        try:
            # Open the input file
            with open(self.InputFile,'r') as infid:
                reader = csv.DictReader(infid, delimiter='\t')
                self.inputNames = reader.fieldnames 
                inData = {}
                for name in self.inputNames:
                    inData[name] = np.array([])
                
                # Check that the input file contains all the needed variables
                success=self.CheckDataPointSeriers()
                if not success:
                    return
                
                # Loop all the lines in the table
                for dataRow in reader:
                    
                    for dataName in self.inputNames:
                        inData[dataName] = np.append(inData[dataName], float(dataRow[dataName]))
                    
                    # Fill in data fields which might not be in the input file
                    if 'SZA' not in self.inputNames:
                        sza, _ = TSEB.met.Get_SunAngles(self.lat, self.lon, self.stdlon, inData['DOY'][-1], inData['Time'][-1]) 
                        addData(inData, 'SZA', sza)
                    if 'SAA' not in self.inputNames:
                        _ , saa = TSEB.met.Get_SunAngles(self.lat, self.lon, self.stdlon, inData['DOY'][-1], inData['Time'][-1])
                        addData(inData, 'SAA', saa)
                    if 'p' not in self.inputNames: # Estimate barometric pressure from the altitude if not included in the table
                        p = TSEB.met.CalcPressure(self.alt)
                        addData(inData, 'p', p) 
                    if 'fc' not in self.inputNames: # Fractional cover
                        addData(inData, 'fc', self.f_c) # Use default value
                    if 'wc' not in self.inputNames: # Canopy width to height ratio
                        addData(inData, 'wc', self.wc) # Use default value
                    if 'fg' not in self.inputNames: # Green fraction
                        addData(inData, 'fg', self.f_g) # Use default value
                    
                    # Esimate diffuse and direct irradiance
                    difvis, difnir, fvis, fnir = TSEB.rad.CalcDifuseRatio(inData['Sdn'][-1], inData['SZA'][-1], press = inData['p'][-1])
                    addData(inData, 'fvis', fvis)
                    addData(inData, 'fnir', fnir)
                    addData(inData, 'Skyl', difvis*fvis+difnir*fnir)
                    addData(inData, 'Sdn_dir', inData['Sdn'][-1]*(1.0-inData['Skyl'][-1]))
                    addData(inData, 'Sdn_dif', inData['Sdn'][-1]*inData['Skyl'][-1])

                    # Incoming long wave radiation
                    if 'Ldn' not in self.inputNames:
                        # Calculate downwelling LW radiance otherwise
                        emisAtm = TSEB.rad.CalcEmiss_atm(inData['ea'][-1], inData['Ta'][-1])
                        Lsky = emisAtm * TSEB.met.CalcStephanBoltzmann(inData['Ta'][-1])
                        addData(inData, 'Lsky', Lsky)                        
                    else:
                        inData['Lsky'] = inData['Ldn']
                        
                    # Calculate Roughness
                    z_0M, d_0 = TSEB.res.CalcRoughness (inData['LAI'][-1:], inData['hc'][-1:], inData['wc'][-1:], self.LANDCOVER)
                    addData(inData, 'z_0M', z_0M)
                    addData(inData, 'd_0', d_0)
        
        except IOError:
            print('Error reading input file : '+self.InputFile)
            return                        
         
        # get the Soil Heat flux if CalcG includes the option of measured G
        if self.CalcG[0]==0: # Constant G
            if 'G' in self.inputNames:
                self.CalcG[1] = inData['G']
        elif self.CalcG[0] == 2: # Santanello and Friedls G
            self.CalcG[1][0] = inData['Time'] # Set the time in the CalcG flag to compute the Santanello and Friedl G
            self.CalcG[1][1] = np.ones(len(inData['Time']))*self.CalcG[1][1]
            self.CalcG[1][2] = np.ones(len(inData['Time']))*self.CalcG[1][2]
            self.CalcG[1][3] = np.ones(len(inData['Time']))*self.CalcG[1][3]

        #======================================
        # Run the chosen model
        
        if self.TSEB_MODEL=='DTD':
            [flag, Ts, Tc, T_AC,S_nS, S_nC, L_nS,L_nC, LE_C,H_C,LE_S,H_S,G,R_s,R_x,R_a,u_friction, L,Ri,
                 n_iterations]=TSEB.DTD(inData['Trad_0'], inData['Trad'], inData['VZA'], inData['Ta_0'], inData['Ta'],
                                        inData['u'], inData['ea'], inData['p'], inData['Sdn_dir'], inData['Sdn_dif'], 
                                        inData['fvis'], inData['fnir'], inData['SZA'], inData['Lsky'], inData['LAI'], 
                                        inData['hc'], self.emisVeg, self.emisGrd, self.spectraVeg, self.spectraGrd, inData['z_0M'],
                                        inData['d_0'], self.zu,self.zt, f_c=inData['fc'], wc=inData['wc'], f_g=inData['fg'],
                                        leaf_width=self.leaf_width, z0_soil=self.z0_soil, alpha_PT=self.Max_alpha_PT, CalcG=self.CalcG)
        
        elif self.TSEB_MODEL=='TSEB_PT':        
            [flag, Ts, Tc, T_AC,S_nS, S_nC, L_nS,L_nC, LE_C,H_C,LE_S,H_S,G,R_s,R_x,R_a,u_friction, L,
                 n_iterations]=TSEB.TSEB_PT(inData['Trad'], inData['VZA'], inData['Ta'], inData['u'], inData['ea'], inData['p'],
                                            inData['Sdn_dir'], inData['Sdn_dif'], inData['fvis'], inData['fnir'], inData['SZA'],
                                            inData['Lsky'], inData['LAI'], inData['hc'], self.emisVeg, self.emisGrd, self.spectraVeg, 
                                            self.spectraGrd, inData['z_0M'], inData['d_0'], self.zu,self.zt, f_c=inData['fc'], 
                                            f_g=inData['fg'], wc=inData['wc'], leaf_width=self.leaf_width, z0_soil=self.z0_soil,
                                            alpha_PT=self.Max_alpha_PT, CalcG=self.CalcG)
        
        elif self.TSEB_MODEL=='TSEB_2T':
            # Run TSEB with the component temperatures Ts and Tc    
            [flag, T_AC,S_nS, S_nC, L_nS,L_nC, LE_C,H_C,LE_S,H_S,G,R_s,R_x,R_a,u_friction, L,
                 n_iterations] = TSEB.TSEB_2T(
                            inData['Tc'], inData['Ts'], inData['Ta'], inData['u'], inData['ea'], inData['p'],
                            inData['Sdn_dir'], inData['Sdn_dif'], inData['fvis'], inData['fnir'], inData['SZA'],
                            inData['Lsky'], inData['LAI'], inData['hc'], self.emisVeg, self.emisGrd, self.spectraVeg, 
                            self.spectraGrd, inData['z_0M'], inData['d_0'], self.zu, self.zt, f_c=inData['fc'], 
                            f_g=inData['fg'], wc=inData['wc'], leaf_width=self.leaf_width, z0_soil=self.z0_soil,
                            alpha_PT=self.Max_alpha_PT, CalcG=self.CalcG)
            Ts = inData['Ts']
            Tc = inData['Tc'] 
            
        elif self.TSEB_MODEL=='OSEB':  
            if np.nanmean(inData['LAI']) > 0: # Assumption of emissivity
                emis = self.emisVeg
            else:
                emis = self.emisGrd
            # Get kB value
            import src.resistances as res_kB
            if isinstance(self.kB_method, float):
                kB_value = self.kB_method
            elif isinstance(self.kB_method, int):
                kB_value = self.kB_method
            elif isinstance(self.kB_method, str):
                if self.kB_method == 'Lhomme':
                    kB_value = res_kB.Calc_kB_Lhomme(inData['LAI'])
                elif self.kB_method == 'Kustas':
                    kB_value = res_kB.Calc_kB_Kustas(inData['Trad'],
                                                     inData['Ta'],inData['u'])
            else:
                print('Could not calculate kB for OSEB model.')
            albedo=fvis*self.spectraGrd['rsoilvis']+fnir* self.spectraGrd['rsoilnir']
            [flag,S_n, L_n, LE,H,G,R_a,u_friction, L,n_iterations] = TSEB.OSEB(inData['Trad'],
                inData['Ta'],inData['u'],inData['ea'],inData['p'],
                inData['Sdn_dir'] + inData['Sdn_dif'],inData['Lsky'],emis,albedo,inData['z_0M'],inData['d_0'],
                self.zu,self.zt,self.CalcG, kB=kB_value )#, T0_K = [], kB = 0.0)
            
        
        if self.TSEB_MODEL=='OSEB':
            Rn=S_n + L_n
            
            #======================================
            # Save output file        
            
            # Output Headers
            outputTxtFieldNames = ['Year', 'DOY', 'Time','LAI','f_g', 'skyl', 'VZA', 
                                   'SZA', 'SAA','L_sky','Rn_model', 
                                   'LE_model', 'H_model',  
                                   'flag', 'zo', 'd', 'G_model', 'R_a', 
                                   'u_friction', 'L',  'n_iterations']
            
            # Create the ouput directory if it doesn't exist
            outdir=dirname(self.OutputFile)
            if not exists(outdir):
                mkdir(outdir)
            
            # Open output file and write the data
            with open (self.OutputFile, 'w') as fp:
                writer = csv.writer(fp, delimiter='\t')
                writer.writerow(outputTxtFieldNames)
                for row in range(LE.size):
                    outData = [ inData['Year'][row], inData['DOY'][row], inData['Time'][row], inData['LAI'][row], 
                                inData['fg'][row], inData['Skyl'][row], inData['VZA'][row], inData['SZA'][row], 
                                inData['SAA'][row], inData['Lsky'][row], Rn[row],LE[row], H[row], 
                                flag[row], inData['z_0M'][row], inData['d_0'][row], 
                                G[row], R_a[row], u_friction[row], L[row], n_iterations]
                    writer.writerow(outData)
            print('Done')    
        else:
            # Calculate the bulk fluxes
            LE=LE_C+LE_S
            H=H_C+H_S
            Rn=S_nC+S_nS+L_nC+L_nS
            
            
            #======================================
            # Save output file        
            
            # Output Headers
            outputTxtFieldNames = ['Year', 'DOY', 'Time','LAI','f_g', 'skyl', 'VZA', 
                                   'SZA', 'SAA','L_sky','Rn_model','Rn_sw_veg', 'Rn_sw_soil', 
                                   'Rn_lw_veg', 'Rn_lw_soil', 'Tc', 'Ts', 'Tac', 
                                   'LE_model', 'H_model', 'LE_c', 'H_c', 'LE_s', 'H_s', 
                                   'flag', 'zo', 'd', 'G_model', 'R_s', 'R_x', 'R_a', 
                                   'u_friction', 'L',  'n_iterations']
            
            # Create the ouput directory if it doesn't exist
            outdir=dirname(self.OutputFile)
            if not exists(outdir):
                mkdir(outdir)
            
            # Open output file and write the data
            with open (self.OutputFile, 'w') as fp:
                writer = csv.writer(fp, delimiter='\t')
                writer.writerow(outputTxtFieldNames)
                for row in range(LE.size):
                    outData = [ inData['Year'][row], inData['DOY'][row], inData['Time'][row], inData['LAI'][row], 
                                inData['fg'][row], inData['Skyl'][row], inData['VZA'][row], inData['SZA'][row], 
                                inData['SAA'][row], inData['Lsky'][row], Rn[row], S_nC[row], S_nS[row], L_nC[row], 
                                L_nS[row], Tc[row], Ts[row], T_AC[row], LE[row], H[row], LE_C[row], H_C[row], 
                                LE_S[row], H_S[row], flag[row], inData['z_0M'][row], inData['d_0'][row], 
                                G[row], R_s[row], R_x[row], R_a[row], u_friction[row], L[row], n_iterations]
                    writer.writerow(outData)
            print('Done')

