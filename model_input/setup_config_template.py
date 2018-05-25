# -*- coding: utf-8 -*-
"""
Sets up the input data to process TIR data from Petit Nobressart using TSEB, DTD or OSEB 
"""
CONFIG_FILE = {
        
    'ModelMode'     : 'TSEB_PT',       # 'TSEB_PT' or 'OSEB'
    
    'SiteProps' : {
        'lat'       : ,     # Latitude
        'lon'       : ,     # Longitude
        'stdlon'    : ,     # Standard longitude, central longitude of the 
                            # Time zone of the site (degrees)
        'altitude'  : ,     # (m) Altitude
        'zt'        : ,     # (m) Measurement height temperature
        'zu'        : ,     # (m), Measurement height wind
        'landCover' : ,     # Primary land cover CROP=11, GRASS=2, SHRUB=5, 
                            # CONIFER=4, BROADLEAVED=3
        },
        
    'Meteo'         : {     # Meteo will be overwritten by the data in the
        'DOY'       : 0,    # meteorological input file
        'Time'      : 0,
        'Ta_0'      : 0,
        'Ta_1'      : 0,
        'u'         : 0,
        'p'         : 0,
        'ea'        : 0,
        'Sdn'       : 0,
        'Ldn'       : '',
        },
            
    'SiteVegProps'  : {
                            # LAI, Fc, Hc, Wc, and Fg will be overwritten
        'LAI'       : ,     # Effective Leaf Area Index (m2/m2) 
        'Fc'        : ,     # Vegetation Fractional Cover 
        'Hc'        : ,     # Canopy height (m), 
        'Wc'        : ,     # Canopy height/with ratio (wc/hc) 
        'Fg'        : ,     # Green Fraction 
        'xLAD'      : ,     # Cambpbell 1990 leaf inclination distribution 
                                    # parameter:[x_LAD=1 for spherical LIDF, 
                                    # x_LAD=0 for vertical LIDF 
                                    # x_LAD=float(inf) for horzontal LIDF] 
        'leafWidth' : ,     # leaf effective width (m)
        
        # Emissivities
        'emisVeg'   : ,     # Leaf emissivity
        'emisSoil'  : ,     # Soil emissivity
        
        # Canopy & soil spectral properties
        'rhovis'    : ,    # rho_leaf_vis: visible reflectance 0.07
        'tauvis'    : ,    # tau_leaf_vis: visible transmittance 0.08
        'rhonir'    : ,    # rho_leaf_nir: NIR reflectance 0.32
        'taunir'    : ,    # tau_leaf_nir: NIR transmittance 0.33
        'rsoilvis'  : ,    # rsoilv: visible reflectance 0.15
        'rsoilnir'  : ,    # rsoiln: NIR reflectance 0.25
        },
                                   
           
    'AncProps' : {
        'maxAlphaPT': ,    # Initial value for Priestley Taylor canopy transpiration
        'VZA'       : ,    # View Zenith Angle (degrees) 
        'z0soil'    : ,    # (m) Bare soil roughness length 
        'useMask'   : ,    # Processing Mask (boolean) 
        },
            
    'SHFinfo' : {
        'CalcG'     : ,    # switch to select soil heat flux measurment method
                           # 0: Use a constant G, usually use G_Constant=0 
                           #    to ignore the computation of G
                           # 1: default, estimate G as a ratio of Rn_soil, 
                           #    default G_ratio=0.35
                           # 2: estimate G from Santanello and Friedl 
                           #    with GAmp the maximum ration amplitude, 
                           #    Gphase, the time shift between G and Rn 
                           #    (hours) and Gshape the typical diurnal shape 
                           #    (hours)
        'G_ratio'   : ,    # estimate G as a ratio of Rn_soil
        'G_constant': ,    # estimate G as a constant
        
        # estimate G from Santanello and Friedl
        'GAmp'      : ,   # Gamp the maximum ration amplitude 
        'Gphase'    : ,   # Gphase, the time shift between G and Rn (hours) 
        'Gshape'    :     # Gshape the typical diurnal shape (hours)
        },
        
    'PointTimeseriesInput' : {
        'InputFile' : '',
        'kB'        :     # kB value to correct the aerodynamic resistance
                          # 'Kustas', 'Lhomme', or floar 2.3
        }}
            
