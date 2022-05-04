# -*- coding: utf-8 -*-

#################################################################################
#										
# plot_cross_section.py						
#										
# Contains routines to plot neutrino-nucleon cross section data.
#										
# Created: 16/06/2017 16:48
# Last modified: 16/06/2017 16:48
#										
#################################################################################


from numpy import *
from pylab import *
import numpy as np
from matplotlib import *
import matplotlib as mpl
import os
import sys

from global_defs import myLightBlue
from global_defs import myDarkGreen

def Plot_Cross_Section_Measurements(ax):
    conv_fb_to_cm2 = 1.e-39 # 1 fb = 1E-36 cm^2
    scale_factor = 1.e38


    #################################################################################
    #### Cross section data
    #### All of the data is for cross section per nucleon, even if the experiment
    #### was nu + nucleon scattering
    #################################################################################

    #################################################################################
    # T2K (C)
    # PRD 92, 112003 (2015)
    # https://www.hepdata.net/record/ins1329784

    # Table 1 (NUMU)
    # '$E_\nu$ [GeV]','$E_\nu$ [GeV] LOW','$E_\nu$ [GeV] HIGH','$\sigma(E)$ [$10^{-38}$ cm$^{2}$]','error +','error -'
    # 0.43761,0.0,0.6,0.68491,0.13303,-0.11087
    # 0.6482,0.6,0.7,0.7715,0.2071,-0.18971
    # 0.8091,0.7,1.0,1.08959,0.23881,-0.19218
    # 1.19633,1.0,1.5,0.3689,0.34782,-0.32917
    # 3.37397,1.5,30.0,1.37898,0.43842,-0.40189

    arr_raw = [
    [ 0.43761,0.0,0.6,0.68491,0.13303,-0.11087 ],
    [ 0.6482,0.6,0.7,0.7715,0.2071,-0.18971 ],
    [ 1.19633,1.0,1.5,0.3689,0.34782,-0.32917 ],
    [ 3.37397,1.5,30.0,1.37898,0.43842,-0.40189 ]
    ]

    arr_cs_cc_numu_t2k_c_2015_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2015_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2015_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2015_cs_central = [ x[3]/x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2015_cs_hi = [ x[4]/x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2015_cs_lo = [ -x[5]/x[0] for x in arr_raw ]

    # There is no NUMUBAR data
    #################################################################################


    #################################################################################
    # T2K (C)
    # PRD 87, 092003 (2013)
    # https://arxiv.org/abs/1302.4908

    # Eq. (12)
    # '$E_\nu$ [GeV]','$E_\nu$ [GeV] LOW','$E_\nu$ [GeV] HIGH','$\sigma(E)$ [$10^{-38}$ cm$^{2}$]','error +','error -','sys error +'.'sys error -'

    arr_raw = [
    [ 0.85,0.58,1.25,0.691,0.013,0.013,0.084,0.084 ],
    ]

    arr_cs_cc_numu_t2k_c_2013_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2013_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2013_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2013_cs_central = [ x[3]/x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2013_cs_hi = [ sqrt(x[4]**2.0+x[6]**2.0)/x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_c_2013_cs_lo = [ sqrt(x[5]**2.0+x[7]**2.0)/x[0] for x in arr_raw ]

    # There is no NUMUBAR data
    #################################################################################


    #################################################################################
    # T2K (Fe)
    # PRD 90, 052010 (2014)
    # https://arxiv.org/abs/1407.4256

    # Eq. (2)
    # '$E_\nu$ [GeV]','$E_\nu$ [GeV] LOW','$E_\nu$ [GeV] HIGH','$\sigma(E)$ [$10^{-38}$ cm$^{2}$]','error +','error -','sys error +'.'sys error -'

    arr_raw = [
    [ 1.51, 0.6965594138260591, 2.3129181267919714, 1.444, 0.002, 0.002, 0.189, 0.157 ],
    ]

    arr_cs_cc_numu_t2k_fe_2014_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_fe_2014_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_t2k_fe_2014_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_fe_2014_cs_central = [ x[3]/x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_fe_2014_cs_hi = [ sqrt(x[4]**2.0+x[6]**2.0)/x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_fe_2014_cs_lo = [ sqrt(x[5]**2.0+x[7]**2.0)/x[0] for x in arr_raw ]

    # There is no NUMUBAR data
    #################################################################################


    #################################################################################
    # T2K (CH)
    # PRD 90, 052010 (2014)
    # https://arxiv.org/abs/1407.4256

    # Eq. (2)
    # '$E_\nu$ [GeV]','$E_\nu$ [GeV] LOW','$E_\nu$ [GeV] HIGH','$\sigma(E)$ [$10^{-38}$ cm$^{2}$]','error +','error -','sys error +'.'sys error -'

    arr_raw = [
    [ 1.51, 0.7368421052631577, 2.35672514619883, 1.379, 0.009, 0.009, 0.178, 0.147 ],
    ]

    arr_cs_cc_numu_t2k_ch_2014_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_ch_2014_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_t2k_ch_2014_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_ch_2014_cs_central = [ x[3]/x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_ch_2014_cs_hi = [ sqrt(x[4]**2.0+x[6]**2.0)/x[0] for x in arr_raw ]
    arr_cs_cc_numu_t2k_ch_2014_cs_lo = [ sqrt(x[5]**2.0+x[7]**2.0)/x[0] for x in arr_raw ]

    # There is no NUMUBAR data
    #################################################################################


    #################################################################################
    # IHEP-ITEP
    # Sov. J. Nucl. Phys. 30 (1979) 528, 1979
    # https://www.hepdata.net/record/ins141742

    # Table 1 (NUMU)
    # 'E [GEV]','E [GEV] LOW','E [GEV] HIGH','SIG/E [FB/GEV]','error +','error -','sys,DUE TO THE NORMALISATION ERROR +','sys,DUE TO THE NORMALISATION ERROR -'
    # 6.5,5.0,8.0,7.7,0.4,-0.4,'7.0%','-7.0%'
    # 10.0,8.0,12.0,7.4,0.4,-0.4,'7.0%','-7.0%'
    # 16.0,12.0,20.0,7.2,0.4,-0.4,'7.0%','-7.0%'
    # 27.5,20.0,35.0,6.8,0.6,-0.6,'7.0%','-7.0%'

    arr_raw = [  
    [ 6.5,5.0,8.0,7.7,0.4,-0.4,'7.0%','-7.0%' ] ,
    [ 10.0,8.0,12.0,7.4,0.4,-0.4,'7.0%','-7.0%' ],
    [ 16.0,12.0,20.0,7.2,0.4,-0.4,'7.0%','-7.0%' ],
    [ 27.5,20.0,35.0,6.8,0.6,-0.6,'7.0%','-7.0%' ]
    ]

    arr_cs_cc_numu_ihep_itep_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_ihep_itep_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_ihep_itep_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_ihep_itep_cs_central = [ x[3]*conv_fb_to_cm2*scale_factor for x in arr_raw ]
    arr_cs_cc_numu_ihep_itep_cs_hi = [ x[4]*conv_fb_to_cm2*scale_factor for x in arr_raw ]
    arr_cs_cc_numu_ihep_itep_cs_lo = [ -x[5]*conv_fb_to_cm2*scale_factor for x in arr_raw ]

    # Table 2 (NUMUBAR)
    # 'E [GEV]','E [GEV] LOW','E [GEV] HIGH','SIG/E [FB/GEV]','error +','error -','sys,DUE TO THE NORMALISATION ERROR +','sys,DUE TO THE NORMALISATION ERROR -'
    # 6.5,5.0,8.0,3.0,0.2,-0.2,'7.0%','-7.0%'
    # 10.0,8.0,12.0,3.1,0.3,-0.3,'7.0%','-7.0%'
    # 16.0,12.0,20.0,3.2,0.4,-0.4,'7.0%','-7.0%'
    # 27.5,20.0,35.0,3.2,0.6,-0.6,'7.0%','-7.0%'

    arr_raw = [  
    [ 6.5,5.0,8.0,3.0,0.2,-0.2,'7.0%','-7.0%' ] ,
    [ 10.0,8.0,12.0,3.1,0.3,-0.3,'7.0%','-7.0%' ],
    [ 16.0,12.0,20.0,3.2,0.4,-0.4,'7.0%','-7.0%' ],
    [ 27.5,20.0,35.0,3.2,0.6,-0.6,'7.0%','-7.0%' ]
    ]

    arr_cs_cc_numubar_ihep_itep_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numubar_ihep_itep_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numubar_ihep_itep_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numubar_ihep_itep_cs_central = [ x[3]*conv_fb_to_cm2*scale_factor for x in arr_raw ]
    arr_cs_cc_numubar_ihep_itep_cs_hi = [ x[4]*conv_fb_to_cm2*scale_factor for x in arr_raw ]
    arr_cs_cc_numubar_ihep_itep_cs_lo = [ -x[5]*conv_fb_to_cm2*scale_factor for x in arr_raw ]
    #################################################################################


    #################################################################################
    # NOMAD
    # Phys. Lett. B660 (2008) 19-25, 2008
    # https://hepdata.net/record/ins767013

    # Table 1 (NUMU)
    # 'E(P=1) [GEV]','E(P=1) [GEV] LOW','E(P=1) [GEV] HIGH','SIG/E(P=1) [10**-38CM**2/GEV]','stat +','stat -','sys +','sys -'
    # 4.6,2.5,6.0,0.786,0.011,-0.011,0.035,-0.035
    # 6.5,6.0,7.0,0.763,0.011,-0.011,0.036,-0.036
    # 7.5,7.0,8.0,0.722,0.009,-0.009,0.035,-0.035
    # 8.5,8.0,9.0,0.701,0.007,-0.007,0.033,-0.033
    # 9.5,9.0,10.0,0.716,0.007,-0.007,0.033,-0.033
    # 10.5,10.0,11.0,0.706,0.005,-0.005,0.026,-0.026
    # 11.5,11.0,12.0,0.705,0.005,-0.005,0.024,-0.024
    # 12.5,12.0,13.0,0.697,0.005,-0.005,0.024,-0.024
    # 13.5,13.0,14.0,0.7,0.005,-0.005,0.024,-0.024
    # 14.5,14.0,15.0,0.698,0.004,-0.004,0.025,-0.025
    # 16.2,15.0,17.5,0.698,0.003,-0.003,0.025,-0.025
    # 18.7,17.5,20.0,0.7,0.003,-0.003,0.025,-0.025
    # 21.2,20.0,22.5,0.699,0.003,-0.003,0.024,-0.024
    # 23.7,22.5,25.0,0.694,0.003,-0.003,0.024,-0.024
    # 26.2,25.0,27.5,0.694,0.003,-0.003,0.025,-0.025
    # 28.7,27.5,30.0,0.694,0.003,-0.003,0.025,-0.025
    # 32.3,30.0,35.0,0.677,0.003,-0.003,0.026,-0.026
    # 37.3,35.0,40.0,0.681,0.003,-0.003,0.026,-0.026
    # 42.4,40.0,45.0,0.675,0.003,-0.003,0.028,-0.028
    # 47.4,45.0,50.0,0.682,0.004,-0.004,0.027,-0.027
    # 54.6,50.0,60.0,0.67,0.003,-0.003,0.028,-0.028
    # 64.7,60.0,70.0,0.675,0.003,-0.003,0.031,-0.031
    # 74.8,70.0,80.0,0.684,0.003,-0.003,0.037,-0.037
    # 84.8,80.0,90.0,0.678,0.004,-0.004,0.041,-0.041
    # 94.8,90.0,100.0,0.677,0.004,-0.004,0.043,-0.043
    # 107.0,100.0,115.0,0.674,0.004,-0.004,0.048,-0.048
    # 122.0,115.0,130.0,0.661,0.005,-0.005,0.048,-0.048
    # 136.9,130.0,145.0,0.671,0.006,-0.006,0.054,-0.054
    # 165.9,145.0,200.0,0.667,0.004,-0.004,0.054,-0.054
    # 228.3,200.0,300.0,0.721,0.008,-0.008,0.06,-0.06

    arr_raw = [  
    [ 4.6,2.5,6.0,0.786,0.011,-0.011,0.035,-0.035 ],
    [ 6.5,6.0,7.0,0.763,0.011,-0.011,0.036,-0.036 ],
    [ 7.5,7.0,8.0,0.722,0.009,-0.009,0.035,-0.035 ],
    [ 8.5,8.0,9.0,0.701,0.007,-0.007,0.033,-0.033 ],
    [ 9.5,9.0,10.0,0.716,0.007,-0.007,0.033,-0.033 ],
    [ 10.5,10.0,11.0,0.706,0.005,-0.005,0.026,-0.026 ],
    [ 11.5,11.0,12.0,0.705,0.005,-0.005,0.024,-0.024 ],
    [ 12.5,12.0,13.0,0.697,0.005,-0.005,0.024,-0.024 ],
    [ 13.5,13.0,14.0,0.7,0.005,-0.005,0.024,-0.024 ],
    [ 14.5,14.0,15.0,0.698,0.004,-0.004,0.025,-0.025 ],
    [ 16.2,15.0,17.5,0.698,0.003,-0.003,0.025,-0.025 ],
    [ 18.7,17.5,20.0,0.7,0.003,-0.003,0.025,-0.025 ],
    [ 21.2,20.0,22.5,0.699,0.003,-0.003,0.024,-0.024 ],
    [ 23.7,22.5,25.0,0.694,0.003,-0.003,0.024,-0.024 ],
    [ 26.2,25.0,27.5,0.694,0.003,-0.003,0.025,-0.025 ],
    [ 28.7,27.5,30.0,0.694,0.003,-0.003,0.025,-0.025 ],
    [ 32.3,30.0,35.0,0.677,0.003,-0.003,0.026,-0.026 ],
    [ 37.3,35.0,40.0,0.681,0.003,-0.003,0.026,-0.026 ],
    [ 42.4,40.0,45.0,0.675,0.003,-0.003,0.028,-0.028 ],
    [ 47.4,45.0,50.0,0.682,0.004,-0.004,0.027,-0.027 ],
    [ 54.6,50.0,60.0,0.67,0.003,-0.003,0.028,-0.028 ],
    [ 64.7,60.0,70.0,0.675,0.003,-0.003,0.031,-0.031 ],
    [ 74.8,70.0,80.0,0.684,0.003,-0.003,0.037,-0.037 ],
    [ 84.8,80.0,90.0,0.678,0.004,-0.004,0.041,-0.041 ],
    [ 94.8,90.0,100.0,0.677,0.004,-0.004,0.043,-0.043 ],
    [ 107.0,100.0,115.0,0.674,0.004,-0.004,0.048,-0.048 ],
    [ 122.0,115.0,130.0,0.661,0.005,-0.005,0.048,-0.048 ],
    [ 136.9,130.0,145.0,0.671,0.006,-0.006,0.054,-0.054 ],
    [ 165.9,145.0,200.0,0.667,0.004,-0.004,0.054,-0.054 ],
    [ 228.3,200.0,300.0,0.721,0.008,-0.008,0.06,-0.06 ]
    ]

    arr_cs_cc_numu_nomad_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_nomad_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_nomad_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_nomad_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_nomad_cs_hi = [ sqrt( x[4]**2.0 + x[6]**2.0 ) for x in arr_raw ]
    arr_cs_cc_numu_nomad_cs_lo = [ sqrt( x[5]**2.0 + x[7]**2.0 ) for x in arr_raw ]

    # There is no NUMUBAR data
    #################################################################################


    #################################################################################
    # MINOS
    # Phys. Rev. D81 (2010) 072002 
    # https://arxiv.org/abs/0910.2201

    # Table III (NUMU)
    # 'E [GEV]','E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -','sys ERROR +','sys ERROR -','sys,DUE TO THE NORMALISATION ERROR +','sys,DUE TO THE NORMALISATION ERROR -', 'total error+', total error-'

    arr_raw = [
    [ 3.48, 3.0, 4.0, 0.748, 0.003, 0.003, 0.058, 0.058, 0.017, 0.017, 0.061, 0.061 ],
    [ 4.45, 4.0, 5.0, 0.711, 0.004, 0.004, 0.029, 0.029, 0.017, 0.017, 0.033, 0.033 ],
    [ 5.89, 5.0, 7.0, 0.708, 0.005, 0.005, 0.027, 0.027, 0.016, 0.016, 0.032, 0.032 ],
    [ 7.97, 7.0, 9.0, 0.722, 0.006, 0.006, 0.041, 0.041, 0.017, 0.017, 0.045, 0.045 ],
    [ 10.45, 9.0, 12.0, 0.699, 0.005, 0.005, 0.041, 0.041, 0.014, 0.014, 0.043, 0.043 ],
    [ 13.43, 12.0, 15.0, 0.691, 0.006, 0.006, 0.023, 0.023, 0.014, 0.014, 0.028, 0.028 ],
    [ 16.42, 15.0, 18.0, 0.708, 0.008, 0.008, 0.012, 0.012, 0.014, 0.014, 0.020, 0.020 ],
    [ 19.87, 18.0, 22.0, 0.689, 0.006, 0.006, 0.009, 0.009, 0.012, 0.012, 0.016, 0.016 ],
    [ 23.88, 22.0, 26.0, 0.683, 0.008, 0.008, 0.005, 0.005, 0.012, 0.012, 0.015, 0.015 ],
    [ 27.89, 26.0, 30.0, 0.686, 0.010, 0.010, 0.004, 0.004, 0.012, 0.012, 0.016, 0.016 ],
    [ 32.81, 30.0, 36.0, 0.675, 0.010, 0.010, 0.002, 0.002, 0.011, 0.011, 0.016, 0.016 ],
    [ 38.87, 36.0, 42.0, 0.675, 0.013, 0.013, 0.005, 0.005, 0.011, 0.011, 0.018, 0.018 ],
    [ 45.77, 42.0, 50.0, 0.676, 0.014, 0.014, 0.004, 0.004, 0.011, 0.011, 0.019, 0.019 ]
    ]

    arr_cs_cc_numu_minos_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_minos_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_minos_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_minos_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_minos_cs_hi = [ x[10] for x in arr_raw ]
    arr_cs_cc_numu_minos_cs_lo = [ x[11] for x in arr_raw ]

    # Table III (NUMUBAR)
    # 'E [GEV]','E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -','sys ERROR +','sys ERROR -','sys,DUE TO THE NORMALISATION ERROR +','sys,DUE TO THE NORMALISATION ERROR -', 'total error+', total error-'

    arr_raw = [
    [ 6.07, 5.0, 7.0, 0.305, 0.005, 0.005, 0.027, 0.027, 0.007, 0.007, 0.029, 0.029 ],
    [ 7.99, 7.0, 9.0, 0.300, 0.005, 0.005, 0.021, 0.021, 0.007, 0.007, 0.022, 0.022 ],
    [ 10.43, 9.0, 12.0, 0.303, 0.004, 0.004, 0.018, 0.018, 0.006, 0.006, 0.019, 0.019 ],
    [ 13.42, 12.0, 15.0, 0.314, 0.005, 0.005, 0.014, 0.014, 0.006, 0.006, 0.016, 0.016 ],
    [ 16.41, 15.0, 18.0, 0.304, 0.007, 0.007, 0.007, 0.007, 0.006, 0.006, 0.012, 0.012 ],
    [ 19.82, 18.0, 22.0, 0.316, 0.006, 0.006, 0.011, 0.011, 0.005, 0.005, 0.013, 0.013 ],
    [ 23.82, 22.0, 26.0, 0.320, 0.009, 0.009, 0.004, 0.004, 0.005, 0.005, 0.011, 0.011 ],
    [ 27.84, 26.0, 30.0, 0.332, 0.012, 0.012, 0.005, 0.005, 0.006, 0.006, 0.015, 0.015 ],
    [ 32.72, 30.0, 36.0, 0.325, 0.014, 0.014, 0.006, 0.006, 0.005, 0.005, 0.016, 0.016 ],
    [ 38.74, 36.0, 42.0, 0.352, 0.021, 0.021, 0.011, 0.011, 0.006, 0.006, 0.024, 0.024 ],
    [ 45.61, 42.0, 50.0, 0.324, 0.023, 0.023, 0.013, 0.013, 0.005, 0.005, 0.027, 0.027 ]
    ]

    arr_cs_cc_numubar_minos_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numubar_minos_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numubar_minos_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numubar_minos_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numubar_minos_cs_hi = [ x[10] for x in arr_raw ]
    arr_cs_cc_numubar_minos_cs_lo = [ x[11] for x in arr_raw ]
    #################################################################################

    #################################################################################
    # SKAT
    # Phys. Lett. B81 (1979) 255 
    # http://www.sciencedirect.com/science/article/pii/0370269379905367?via%3Dihub

    # Table 1 (NUMU)
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 3.0, 5.0, 0.68, 0.07, 0.07 ],
    [ 5.0, 9.0, 0.70, 0.09, 0.09 ],
    [ 9.0, 16.0, 0.86, 0.13, 0.13 ],
    [ 16.0, 30.0, 0.70, 0.15, 0.15 ],
    ]

    arr_cs_cc_numu_skat_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numu_skat_energy_lo = [ x[0]-(x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numu_skat_energy_hi = [ (x[0]+x[1])/2.0-x[1] for x in arr_raw ]
    arr_cs_cc_numu_skat_cs_central = [ x[2] for x in arr_raw ]
    arr_cs_cc_numu_skat_cs_hi = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_skat_cs_lo = [ x[4] for x in arr_raw ]

    # There is no NUMUBAR data
    #################################################################################

    #################################################################################
    # SciBooNE
    # Phys. Rev. D83 (2011) 012005 
    # https://arxiv.org/abs/1011.2131

    # Table XII (NUMU) -- NUANCE-based
    # 'E [GEV]','E [GEV] LOW','E [GEV] HIGH','SIG [10**-38CM**2]','error +','error -'

    arr_raw = [
    [ 0.38, 0.25, 0.50, 0.340, 0.096, 0.096 ],
    [ 0.62, 0.50, 0.75, 0.639, 0.081, 0.081 ],
    [ 0.87, 0.75, 1.00, 1.01, 0.09, 0.09 ],
    [ 1.11, 1.00, 1.25, 1.29, 0.15, 0.15 ],
    [ 1.43, 1.25, 1.75, 1.56, 0.28, 0.28 ],
    [ 2.47, 1.75, 3.00, 1.66, 0.37, 0.37 ]
    ]

    arr_cs_cc_numu_sciboone_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_sciboone_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_sciboone_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_sciboone_cs_central = [ x[3]/x[0] for x in arr_raw ]
    arr_cs_cc_numu_sciboone_cs_hi = [ x[4]/x[0] for x in arr_raw ]
    arr_cs_cc_numu_sciboone_cs_lo = [ x[5]/x[0] for x in arr_raw ]

    # There is no NUMUBAR data
    #################################################################################

    #################################################################################
    # GGM-PS
    # Phys. Lett. B84 (1979) 281-284, 1979
    # https://hepdata.net/record/ins141175

    # Table 1 (NUMU)
    # 'PLAB [GEV]','PLAB [GEV] LOW','PLAB [GEV] HIGH','SIG/E [10**-38 CM**2/GEV]','error +','error -'
    # 2.87,1.5,15.0,0.69,0.05,-0.05
    # 9.05,6.0,15.0,0.61,0.06,-0.06

    arr_raw = [
    [ 2.87,1.5,15.0,0.69,0.05,-0.05 ],
    [ 9.05,6.0,15.0,0.61,0.06,-0.06 ]
    ]

    arr_cs_cc_numu_ggm_ps_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_ggm_ps_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_ggm_ps_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_ggm_ps_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_ggm_ps_cs_hi = [ x[4] for x in arr_raw ]
    arr_cs_cc_numu_ggm_ps_cs_lo = [ -x[5] for x in arr_raw ]

    # Table 2 (NUMUBAR)
    # 'PLAB [GEV]','PLAB [GEV] LOW','PLAB [GEV] HIGH','SIG/E [10**-38 CM**2/GEV]','error +','error -'
    # 3.0,1.5,15.0,0.26,0.02,-0.02

    arr_raw = [
    [ 3.0,1.5,15.0,0.26,0.02,-0.02 ]
    ]

    arr_cs_cc_numubar_ggm_ps_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_ps_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_ps_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_ps_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_ps_cs_hi = [ x[4] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_ps_cs_lo = [ -x[5] for x in arr_raw ]
    #################################################################################

    #################################################################################
    # ANL
    # Phys. Rev. D19 (1979) 2521, 1979
    # https://hepdata.net/record/ins7237

    # Table 1 (NUMU)
    # 'PLAB [GEV]','PLAB [GEV] LOW','PLAB [GEV] HIGH','SIG/E [10**-38 CM**2/GEV]','error +','error -'
    # 3.1,0.2,6.0,0.87,0.03,-0.03
    # 1.1,1.1,1.1,0.76,0.06,-0.06

    arr_raw = [
    [ 3.1,0.2,6.0,0.87,0.03,-0.03 ],
    [ 1.1,1.1,1.1,0.76,0.06,-0.06 ]
    ]

    arr_cs_cc_numu_anl_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_anl_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_anl_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_anl_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_anl_cs_hi = [ x[4] for x in arr_raw ]
    arr_cs_cc_numu_anl_cs_lo = [ -x[5] for x in arr_raw ]

    # There is no NUMUBAR data
    #################################################################################

    #################################################################################
    # BEBC
    # Z. Phys. C1 (1979) 143 
    # http://dx.doi.org/10.1007/BF01445406

    # Table 2 (NUMU)
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 30.0, 90.0, 0.62, 0.04, 0.04 ],
    [ 90.0, 190.0, 0.63, 0.05, 0.05 ]
    ]

    arr_cs_cc_numu_bebc_1_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numu_bebc_1_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    arr_cs_cc_numu_bebc_1_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numu_bebc_1_cs_central = [ x[2] for x in arr_raw ]
    arr_cs_cc_numu_bebc_1_cs_hi = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_bebc_1_cs_lo = [ x[4] for x in arr_raw ]

    # Table 2 (NUMUBAR)
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 30.0, 90.0, 0.30, 0.02, 0.02 ],
    [ 90.0, 190.0, 0.31, 0.04, 0.04 ]
    ]

    arr_cs_cc_numubar_bebc_1_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numubar_bebc_1_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    arr_cs_cc_numubar_bebc_1_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numubar_bebc_1_cs_central = [ x[2] for x in arr_raw ]
    arr_cs_cc_numubar_bebc_1_cs_hi = [ x[3] for x in arr_raw ]
    arr_cs_cc_numubar_bebc_1_cs_lo = [ x[4] for x in arr_raw ]
    #################################################################################


    #################################################################################
    # BEBC
    # Z. Phys. C2 (1979) 187
    # http://dx.doi.org/10.1007/BF01474659

    # Taken from Table 1 of Z. Phys. C 70 (1996) 39 (NUMU)
    # 'E [GEV]','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 19.0, 0.74, 0.11, 0.11 ],
    [ 34.0, 0.73, 0.10, 0.10 ]
    ]

    arr_cs_cc_numu_bebc_2_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_bebc_2_cs_central = [ x[1] for x in arr_raw ]
    arr_cs_cc_numu_bebc_2_cs_hi = [ x[2] for x in arr_raw ]
    arr_cs_cc_numu_bebc_2_cs_lo = [ x[3] for x in arr_raw ]

    # There is no NUMUBAR data
    #################################################################################


    #################################################################################
    # IHEP-JINR
    # Z. Phys. C70 (1996) 39
    # http://dx.doi.org/10.1007/s002880050078

    # Table 5 (NUMU)
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -','sys error +','sys error -'

    arr_raw = [
    [ 3.0, 5.0, 0.792, 0.016, 0.016, 0.039, 0.039 ],
    [ 5.0, 7.0, 0.817, 0.013, 0.013, 0.035, 0.035 ],
    [ 7.0, 9.0, 0.779, 0.014, 0.014, 0.035, 0.035 ],
    [ 9.0, 11.0, 0.738, 0.016, 0.016, 0.026, 0.026 ],
    [ 11.0, 13.0, 0.717, 0.018, 0.018, 0.027, 0.027 ],
    [ 13.0, 17.0, 0.683, 0.015, 0.015, 0.029, 0.029 ],
    [ 17.0, 21.0, 0.654, 0.020, 0.020, 0.030, 0.030 ],
    [ 21.0, 25.0, 0.635, 0.024, 0.024, 0.041, 0.041 ],
    [ 25.0, 30.0, 0.609, 0.031, 0.031, 0.054, 0.054 ]
    ]

    arr_cs_cc_numu_ihep_jinr_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numu_ihep_jinr_energy_lo = [ x[0]-(x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numu_ihep_jinr_energy_hi = [ (x[0]+x[1])/2.0-x[1] for x in arr_raw ]
    arr_cs_cc_numu_ihep_jinr_cs_central = [ x[2] for x in arr_raw ]
    arr_cs_cc_numu_ihep_jinr_cs_hi = [ sqrt(x[3]**2.0+x[5]**2.0) for x in arr_raw ]
    arr_cs_cc_numu_ihep_jinr_cs_lo = [ sqrt(x[4]**2.0+x[6]**2.0) for x in arr_raw ]

    # Table 5 (NUMUBAR)
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -','sys error +','sys error -'

    arr_raw = [
    [ 3.0, 5.0, 0.345, 0.016, 0.016, 0.023, 0.023 ],
    [ 5.0, 7.0, 0.347, 0.013, 0.013, 0.020, 0.020 ],
    [ 7.0, 10.0, 0.336, 0.012, 0.012, 0.019, 0.019 ],
    [ 10.0, 13.0, 0.346, 0.015, 0.015, 0.022, 0.022 ],
    [ 13.0, 17.0, 0.323, 0.017, 0.017, 0.022, 0.022 ],
    [ 17.0, 21.0, 0.284, 0.024, 0.024, 0.021, 0.021 ],
    [ 21.0, 26.0, 0.269, 0.029, 0.029, 0.023, 0.023 ],
    [ 26.0, 30.0, 0.279, 0.039, 0.039, 0.031, 0.031 ]
    ]

    arr_cs_cc_numubar_ihep_jinr_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numubar_ihep_jinr_energy_lo = [ x[0]-(x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_numubar_ihep_jinr_energy_hi = [ (x[0]+x[1])/2.0-x[1] for x in arr_raw ]
    arr_cs_cc_numubar_ihep_jinr_cs_central = [ x[2] for x in arr_raw ]
    arr_cs_cc_numubar_ihep_jinr_cs_hi = [ sqrt(x[3]**2.0+x[5]**2.0) for x in arr_raw ]
    arr_cs_cc_numubar_ihep_jinr_cs_lo = [ sqrt(x[4]**2.0+x[6]**2.0) for x in arr_raw ]
    #################################################################################


    #################################################################################
    # CCFR
    # Seligman thesis 1997

    # NUMU
    # 'E [GEV]','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 76.06192263888305, 0.6993630573248406, 0.0, 0.0 ],
    [ 85.99780655502596, 0.7006369426751591, 0.0, 0.0 ],
    [ 95.92525414434554, 0.6968152866242037, 0.0, 0.0 ],
    [ 110.18897372084193, 0.7121019108280253, 0.0, 0.0 ],
    [ 130.04597798118698, 0.705732484076433, 0.0, 0.0 ],
    [ 150.91323237862235, 0.7095541401273884, 0.0, 0.0 ],
    [ 170.76812755726155, 0.7019108280254776, 0.0, 0.0 ],
    [ 190.93727591006876, 0.684076433121019, 0.0, 0.0 ],
    [ 216.079638925212, 0.6700636942675158, 0.0, 0.0 ],
    [ 244.89813135360868, 0.6764331210191081, 0.0, 0.0 ],
    [ 275.024254439617, 0.6726114649681527, 0.0, 0.0 ],
    [ 303.8469650314253, 0.681528662420382, 0.6993630573248406, 0.6649681528662419 ],
    [ 338.9484118614755, 0.6828025477707005, 0.7019108280254776, 0.6649681528662419 ]
    ]

    arr_cs_cc_numu_ccfr_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_ccfr_cs_central = [ x[1] for x in arr_raw ]
    arr_cs_cc_numu_ccfr_cs_hi = [ x[2]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_ccfr_cs_lo = [ x[1]-x[3] for x in arr_raw ]

    # NUMUBAR
    # 'E [GEV]','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 35.73952628138338, 0.34279666562548206, 0.0, 0.0 ],
    [ 45.65740076770575, 0.3350318471337579, 0.0, 0.0 ],
    [ 55.59539376555449, 0.3375796178343947, 0.0, 0.0 ],
    [ 65.52495043657993, 0.3350318471337579, 0.0, 0.0 ],
    [ 75.45872527101702, 0.3350318471337579, 0.0, 0.0 ],
    [ 85.39039102374828, 0.3337579617834393, 0.0, 0.0 ],
    [ 95.32627493989118, 0.3350318471337579, 0.0, 0.0 ],
    [ 109.56890369932931, 0.3375796178343947, 0.0, 0.0 ],
    [ 129.4617623486734, 0.3528662420382165, 0.0, 0.0 ],
    [ 150.31003501075634, 0.34522292993630566, 0.0, 0.0 ],
    [ 169.5068966971781, 0.34012738853503177, 0.0, 0.0 ],
    [ 190.36360568608427, 0.3375796178343947, 0.0, 0.0 ],
    [ 214.52946387143038, 0.3337579617834393, 0.0, 0.0 ],
    [ 243.3289745644746, 0.3286624203821654, 0.0, 0.0 ],
    [ 273.4303998332751, 0.3484849539767412, 0.3659053370524725, 0.3314200889229637 ],
    [ 304.42521440287425, 0.31506625991309356, 0.3371083772742228, 0.29337966057391796 ],
    [ 334.5388857955451, 0.36341671089879657, 0.3996795491380738, 0.32786490870342666 ]
    ]

    arr_cs_cc_numubar_ccfr_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numubar_ccfr_cs_central = [ x[1] for x in arr_raw ]
    arr_cs_cc_numubar_ccfr_cs_hi = [ x[2]-x[1] for x in arr_raw ]
    arr_cs_cc_numubar_ccfr_cs_lo = [ x[1]-x[3] for x in arr_raw ]
    #################################################################################


    #################################################################################
    # CDHS
    # Z. Phys. C 35 (1987) 443 

    # NUMU
    # http://hepdata.cedar.ac.uk/view/ins246156/d1/plain.txt;jsessionid=sm1j32urlv2j
    # Path: /HepData/6594/d1-x1-y1
    # Measured charged current total cross section
    # RE : NUMU NUCLEON --> MU- X
    # SQRT(S) : 7.683 - 16.406 GeV
    # x : PLAB IN GEV
    # y : SIG/E IN 10**-38 CM**2/GEV
    # xdesc	x	xlow	xhigh	y	dy+	dy-	dy+	dy-	
    # 	50.0	10.0	100.0	0.691	+0.007	-0.007	+0.027	-0.027
    # 	85.0	20.0	160.0	0.707	+0.002	-0.002	+0.022	-0.022
    # 	112.0	20.0	200.0	0.708	+0.007	-0.007	+0.029	-0.029
    # 	31.0	31.0	31.0	0.682	+0.008	-0.008	+0.02	-0.02
    # 	50.0	50.0	50.0	0.706	+0.002	-0.002	+0.018	-0.018
    # 	61.0	61.0	61.0	0.707	+0.009	-0.009	+0.023	-0.023
    # 	83.0	83.0	83.0	0.724	+0.012	-0.012	+0.051	-0.051
    # 	121.0	121.0	121.0	0.708	+0.003	-0.003	+0.023	-0.023
    #   143.0	143.0	143.0	0.711	+0.01	-0.01	+0.04	-0.04

    arr_raw = [
    [ 50.0, 10.0,	100.0,	0.691,	+0.007,	-0.007,	+0.027,	-0.027 ],
    [ 85.0,	20.0,	160.0,	0.707,	+0.002,	-0.002,	+0.022,	-0.022 ],
    [ 112.0,20.0,	200.0,	0.708,	+0.007,	-0.007,	+0.029,	-0.029 ],
    [ 31.0,	31.0,	31.0,	0.682,	+0.008,	-0.008,	+0.02,	-0.02 ],
    [ 50.0,	50.0,	50.0,	0.706,	+0.002,	-0.002,	+0.018,	-0.018 ],
    [ 61.0,	61.0,	61.0,	0.707,	+0.009,	-0.009,	+0.023,	-0.023 ],
    [ 83.0,	83.0,	83.0,	0.724,	+0.012,	-0.012,	+0.051,	-0.051 ],
    [ 121.0,121.0,	121.0,	0.708,	+0.003,	-0.003,	+0.023,	-0.023 ],
    [ 143.0,143.0,	143.0,	0.711,	+0.01,	-0.01,	+0.04,  -0.04 ]
    ]

    arr_cs_cc_numu_cdhs_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_cdhs_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_cdhs_energy_hi = [ x[2]-x[0]for x in arr_raw ]
    arr_cs_cc_numu_cdhs_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_cdhs_cs_hi = [ sqrt(x[4]**2.0+x[6]**2.0) for x in arr_raw ]
    arr_cs_cc_numu_cdhs_cs_lo = [ sqrt(x[5]**2.0+x[7]**2.0) for x in arr_raw ]

    # NUMUBAR
    # http://hepdata.cedar.ac.uk/view/ins246156/d2/plain.txt;jsessionid=sm1j32urlv2j
    # Path: /HepData/6594/d2-x1-y1
    # Measured charged current total cross section
    # RE : NUMUBAR NUCLEON --> MU+ X
    # SQRT(S) : 7.683 - 16.406 GeV
    # x : PLAB IN GEV
    # y : SIG/E IN 10**-38 CM**2/GEV
    # xdesc	x	xlow	xhigh	y	dy+	dy-	dy+	dy-	
    # 	47.0	10.0	100.0	0.332	+0.004	-0.004	+0.012	-0.012
    # 	72.0	20.0	160.0	0.333	+0.004	-0.004	+0.009	-0.009
    # 	87.0	20.0	200.0	0.325	+0.007	-0.007	+0.012	-0.012
    # 	31.0	31.0	31.0	0.327	+0.004	-0.004	+0.008	-0.008
    # 	50.0	50.0	50.0	0.332	+0.005	-0.005	+0.007	-0.007
    # 	61.0	61.0	61.0	0.338	+0.009	-0.009	+0.01	-0.01
    # 	83.0	83.0	83.0	0.346	+0.007	-0.007	+0.022	-0.022
    # 	121.0	121.0	121.0	0.337	+0.008	-0.008	+0.012	-0.012
    # 	143.0	143.0	143.0	0.306	+0.015	-0.015	+0.02	-0.02

    arr_raw = [
    [ 47.0, 10.0,	100.0,	0.332,	+0.004,	-0.004,	+0.012,	-0.012 ],
    [ 72.0, 20.0,	160.0,	0.333,	+0.004,	-0.004,	+0.009,	-0.009 ],
    [ 87.0, 20.0,	200.0,	0.325,	+0.007,	-0.007,	+0.012,	-0.012 ],
    [ 31.0, 31.0,	31.0,	0.327,	+0.004,	-0.004,	+0.008,	-0.008 ],
    [ 50.0, 50.0,	50.0,	0.332,	+0.005,	-0.005,	+0.007,	-0.007 ],
    [ 61.0, 61.0,	61.0,	0.338,	+0.009,	-0.009,	+0.01,	-0.01 ],
    [ 83.0, 83.0,	83.0,	0.346,	+0.007,	-0.007,	+0.022,	-0.022 ],
    [ 121.0,121.0,	121.0,	0.337,	+0.008,	-0.008,	+0.012,	-0.012 ],
    [ 143.0,143.0,	143.0,	0.306,	+0.015,	-0.015,	+0.02,	-0.02 ]
    ]

    arr_cs_cc_numubar_cdhs_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numubar_cdhs_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numubar_cdhs_energy_hi = [ x[2]-x[0]for x in arr_raw ]
    arr_cs_cc_numubar_cdhs_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numubar_cdhs_cs_hi = [ sqrt(x[4]**2.0+x[6]**2.0) for x in arr_raw ]
    arr_cs_cc_numubar_cdhs_cs_lo = [ sqrt(x[5]**2.0+x[7]**2.0) for x in arr_raw ]
    #################################################################################


    #################################################################################
    # BNL
    # Phys. Rev. D 25 (1982) 617-623 

    # NUMU
    # http://hepdata.cedar.ac.uk/view/ins177607
    # *dataset:
    # *dscomment: Measured charged current total cross section.
    # Axis error includes +- 0.0/0.0 contribution (?////SYSTEMATIC ERROR NOT GIVENNEUTRAL CURRENT AND NEUTRAL PARTICLES INDUCED REACTIONS, RESCATTERING IN DEUTERIUM).
    # *reackey: NUMU NUCLEON --> MU- X
    # *obskey: SIG
    # *qual: RE : NUMU NUCLEON --> MU- X
    # *qual: SQRT(S) IN GEV : 1.481 TO 3.867
    # *yheader: SIG/E IN 10**-38 CM**2/GEV
    # *xheader: PLAB IN GEV 
    # *data: x : y 
    #  0.7; 1.09 +- 0.08;
    #  0.9; 0.97 +- 0.07;
    #  1.2; 0.98 +- 0.07;
    #  1.3; 0.85 +- 0.07;
    #  1.6; 0.78 +- 0.06;
    #  1.8; 0.79 +- 0.07;
    #  2.0; 0.8 +- 0.08;
    #  2.7; 0.83 +- 0.08;
    #  3.0; 0.79 +- 0.1;
    #  3.5; 0.91 +- 0.14;
    #  4.0; 0.85 +- 0.17;
    #  5.0; 0.78 +- 0.11;
    #  7.5; 0.7 +- 0.15;
    # *dataend:

    arr_raw = [
    [ 0.7, 1.09, 0.08, 0.08 ],
    [ 0.9, 0.97, 0.07, 0.07 ],
    [ 1.2, 0.98, 0.07, 0.07 ],
    [ 1.3, 0.85, 0.07, 0.07 ],
    [ 1.6, 0.78, 0.06, 0.06 ],
    [ 1.8, 0.79, 0.07, 0.07 ],
    [ 2.0, 0.8, 0.08, 0.08 ],
    [ 2.7, 0.83, 0.08, 0.08 ],
    [ 3.0, 0.79, 0.1, 0.1 ],
    [ 3.5, 0.91, 0.14, 0.14 ],
    [ 4.0, 0.85, 0.17, 0.17 ],
    [ 5.0, 0.78, 0.11, 0.11 ],
    [ 7.5, 0.7, 0.15, 0.15 ]
    ]

    arr_cs_cc_numu_bnl_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_bnl_cs_central = [ x[1] for x in arr_raw ]
    arr_cs_cc_numu_bnl_cs_hi = [ x[2] for x in arr_raw ]
    arr_cs_cc_numu_bnl_cs_lo = [ x[3] for x in arr_raw ]

    # No NUMUBAR data
    #################################################################################


    #################################################################################
    # IceCube 4yr
    # Our work

    # Average between NUMU and NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG [CM**2]','error +','error -'

    # arr_raw = [
    # [ 20.e3, 100.e3, 10.**(-34.21), 10.**(-34.21+0.48), 10.**(-34.21-0.48) ],
    # [ 100.e3, 500.e3, 10.**(-32.92), 10.**(-32.92+0.90), 10.**(-32.92-0.90) ],
    # [ 500.e3, 2000.e3, 10.**(-32.5), 10.**(-32.5+1.2), 10.**(-32.5-1.2) ],
    # ]
    # arr_cs_cc_nuavg_ic_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_central = [ 1.e38 * x[2] / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_hi = [ 1.e38 * (x[3]-x[2]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_lo = [ 1.e38 * (x[2]-x[4]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    #################################################################################


    ################################################################################  
    # IceCube 6yr
    # Our work, eq-nu-nubar-sigma

    # Average between NUMU and NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG [CM**2]','error +','error -'

    # Low-res, obtained on laptop
    # arr_raw = [
    # [ 18.e3, 100.e3, 10.**(-34.17), 10.**(-34.17+0.46), 10.**(-34.17-0.46) ],
    # [ 100.e3, 500.e3, 10.**(-33.46), 10.**(-33.46+0.32), 10.**(-33.46-0.32) ],
    # [ 500.e3, 2000.e3, 10.**(-32.5), 10.**(-32.5+1.2), 10.**(-32.5-1.2) ]
    # ]

    # Hi-res, obtained on ruby
    # arr_raw = [
    # [ 18.e3, 100.e3, 10.**(-34.23), 10.**(-34.23+0.48), 10.**(-34.23-0.48) ],
    # [ 100.e3, 500.e3, 10.**(-33.68), 10.**(-33.68+0.67), 10.**(-33.68-0.67) ],
    # [ 500.e3, 2000.e3, 10.**(-31.8), 10.**(-31.8+1.4), 10.**(-31.8-1.4) ]
    # ]
    # arr_cs_cc_nuavg_ic_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_central = [ 1.e38 * x[2] / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_hi = [ 1.e38 * (x[3]-x[2]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_lo = [ 1.e38 * (x[2]-x[4]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    ################################################################################


    ################################################################################  
    # IceCube 6yr - 3 bins
    # Our work, diff-nu-nubar-sigma

    # Average between NUMU and NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG [CM**2]','error +','error -'

    # Low-res, obtained on laptop
    # arr_raw = [
    # [ 18.e3, 100.e3, 10.**(-34.09), 10.**(-34.09+0.63), 10.**(-34.09-0.63) ],
    # [ 100.e3, 500.e3, 10.**(-33.75), 10.**(-33.75+0.48), 10.**(-33.75-0.48) ],
    # [ 500.e3, 2000.e3, 10.**(-32.43), 10.**(-32.43+1.2), 10.**(-32.43-1.2) ]
    # ]

    # # Hi-res, obtained on ruby
    # # arr_raw = [
    # # ]

    # arr_cs_cc_nuavg_ic_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_central = [ 1.e38 * x[2] / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_hi = [ 1.e38 * (x[3]-x[2]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_lo = [ 1.e38 * (x[2]-x[4]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    ################################################################################


    ################################################################################  
    # IceCube 6yr - 6 bins
    # Our work, diff-nu-nubar-sigma

    # Average between NUMU and NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG [CM**2]','error +','error -'

    # Low-res, obtained on laptop
    # arr_raw = [
    # [ 18.e3, 30.e3, 10.**(-34.15), 10.**(-34.15+0.69), 10.**(-34.15-0.69) ],
    # [ 30.e3, 50.e3, 10.**(-34.23), 10.**(-34.23+0.61), 10.**(-34.23-0.61) ],
    # [ 50.e3, 84.e3, 10.**(-33.69), 10.**(-33.69+0.53), 10.**(-33.69-0.53) ],
    # [ 84.e3, 140.e3, 10.**(-33.20), 10.**(-33.20+0.73), 10.**(-33.20-0.73) ],
    # [ 140.e3, 400.e3, 10.**(-33.76), 10.**(-33.76+0.60), 10.**(-33.76-0.60) ],
    # [ 400.e3, 2004.e3, 10.**(-32.43), 10.**(-32.43+1.20), 10.**(-32.43-1.20) ]
    # ]

    # Hi-res, obtained on ruby
    # arr_raw = [
    # ]

    # arr_cs_cc_nuavg_ic_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_central = [ 1.e38 * x[2] / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_hi = [ 1.e38 * (x[3]-x[2]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_lo = [ 1.e38 * (x[2]-x[4]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    ################################################################################

    ################################################################################  
    # IceCube 6yr - 5 bins
    # Our work, diff-nu-nubar-sigma

    # Average between NUMU and NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG [CM**2]','error +','error -'

    # Low-res, obtained on laptop
    # arr_raw = [
    # [ 18.e3, 50.e3, 10.**(-34.05), 10.**(-34.05+0.60), 10.**(-34.05-0.60) ],
    # [ 50.e3, 100.e3, 10.**(-33.71), 10.**(-33.71+0.47), 10.**(-33.71-0.47) ],
    # [ 100.e3, 200.e3, 10.**(-33.51), 10.**(-33.51+0.46), 10.**(-33.51-0.46) ],
    # [ 200.e3, 400.e3, 10.**(-33.50), 10.**(-33.50+0.80), 10.**(-33.50-0.80) ],
    # [ 400.e3, 2004.e3, 10.**(-32.43), 10.**(-32.43+1.20), 10.**(-32.43-1.20) ]
    # ]

    # Hi-res, obtained on ruby
    # arr_raw = [
    # ]

    # arr_cs_cc_nuavg_ic_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_central = [ 1.e38 * x[2] / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_hi = [ 1.e38 * (x[3]-x[2]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_lo = [ 1.e38 * (x[2]-x[4]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    ################################################################################

    ################################################################################  
    # IceCube 6-yr HESE + 2-yr MESE - 4 bins
    # Our work, diff-nu-nubar-sigma

    # Average between NUMU and NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG [CM**2]','error +','error -'

    # Low-res, obtained on laptop
    # arr_raw = [
    # [ 18.e3, 50.e3, 10.**(-34.05), 10.**(-34.05+0.60), 10.**(-34.05-0.60) ],
    # [ 50.e3, 100.e3, 10.**(-33.71), 10.**(-33.71+0.47), 10.**(-33.71-0.47) ],
    # [ 100.e3, 400.e3, 10.**(-33.60), 10.**(-33.60+0.52), 10.**(-33.60-0.52) ],
    # [ 400.e3, 2004.e3, 10.**(-32.43), 10.**(-32.43+1.20), 10.**(-32.43-1.20) ]
    # ]

    # Hi-res, obtained on ruby
    # arr_raw = [
    # [ 1.e3, 18.e3, 10.**(-34.57), 10.**(-34.57+0.29), 10.**(-34.57-0.29) ], # lo-res for now
    # [ 18.e3, 50.e3, 10.**(-34.32), 10.**(-34.32+0.50), 10.**(-34.32-0.50) ],
    # [ 50.e3, 100.e3, 10.**(-33.90), 10.**(-33.90+0.62), 10.**(-33.90-0.62) ],
    # [ 100.e3, 400.e3, 10.**(-33.71), 10.**(-33.71+0.72), 10.**(-33.71-0.72) ],
    # [ 400.e3, 2004.e3, 10.**(-31.83), 10.**(-31.83+1.40), 10.**(-31.83-1.40) ]
    # ]

    # arr_cs_cc_nuavg_ic_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_central = [ 1.e38 * x[2] / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_hi = [ 1.e38 * (x[3]-x[2]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    # arr_cs_cc_nuavg_ic_cs_lo = [ 1.e38 * (x[2]-x[4]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    ################################################################################

    ################################################################################  
    # IceCube 6-yr HESE - 4 bins -- including self-veto
    # Our work, diff-nu-nubar-sigma

    # Average between NUMU and NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG [CM**2]','error +','error -'

    # Hi-res, obtained on ruby
    arr_raw = [
    [ 18.e3, 50.e3, 10.**(-34.35), 10.**(-34.35+0.53), 10.**(-34.35-0.53) ],
    [ 50.e3, 100.e3, 10.**(-33.80), 10.**(-33.80+0.67), 10.**(-33.80-0.67) ],
    [ 100.e3, 400.e3, 10.**(-33.84), 10.**(-33.84+0.67), 10.**(-33.84-0.67) ],
    [ 400.e3, 2004.e3, 10.**(-31.71), 10.**(-31.71+1.50), 10.**(-31.71-1.50) ]
    ]

    arr_cs_cc_nuavg_ic_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_nuavg_ic_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    arr_cs_cc_nuavg_ic_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_nuavg_ic_cs_central = [ 1.e38 * x[2] / ((x[0]+x[1])/2.0) for x in arr_raw ]
    arr_cs_cc_nuavg_ic_cs_hi = [ 1.e38 * (x[3]-x[2]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    arr_cs_cc_nuavg_ic_cs_lo = [ 1.e38 * (x[2]-x[4]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    ################################################################################
    
    ################################################################################  
    # IceCube 6-yr HESE - 4 bins -- including self-veto
    # Our work, diff-nu-nubar-sigma

    # Average between NUMU and NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG [CM**2]','error +','error -'

    # Hi-res, obtained on ruby
    arr_raw = [
    [ 18.e3, 50.e3, 10.**(-34.35), 10.**(-34.35+0.53), 10.**(-34.35-0.53) ],
    [ 50.e3, 100.e3, 10.**(-33.80), 10.**(-33.80+0.67), 10.**(-33.80-0.67) ],
    [ 100.e3, 400.e3, 10.**(-33.84), 10.**(-33.84+0.67), 10.**(-33.84-0.67) ],
    [ 400.e3, 2004.e3, 10.**(-31.71), 10.**(-31.71+1.50), 10.**(-31.71-1.50) ]
    ]

    arr_cs_cc_nuavg_ic_energy_central = [ (x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_nuavg_ic_energy_lo = [ (x[0]+x[1])/2.0-x[0] for x in arr_raw ]
    arr_cs_cc_nuavg_ic_energy_hi = [ x[1]-(x[0]+x[1])/2.0 for x in arr_raw ]
    arr_cs_cc_nuavg_ic_cs_central = [ 1.e38 * x[2] / ((x[0]+x[1])/2.0) for x in arr_raw ]
    arr_cs_cc_nuavg_ic_cs_hi = [ 1.e38 * (x[3]-x[2]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    arr_cs_cc_nuavg_ic_cs_lo = [ 1.e38 * (x[2]-x[4]) / ((x[0]+x[1])/2.0) for x in arr_raw ]
    ################################################################################



    #################################################################################
    # ArgoNeuT 2014
    # Phys. Rev. D 89, 112003 (2014)
    # https://arxiv.org/abs/1404.4809

    # NUMU
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 9.6,6.5,6.5,0.66,0.03,0.03,0.08,0.08 ]
    ]

    arr_cs_cc_numu_argoneut_14_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_argoneut_14_energy_lo = [ x[1] for x in arr_raw ]
    arr_cs_cc_numu_argoneut_14_energy_hi = [ x[2] for x in arr_raw ]
    arr_cs_cc_numu_argoneut_14_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_argoneut_14_cs_hi = [ sqrt(x[4]**2.0+x[6]**2.0) for x in arr_raw ]
    arr_cs_cc_numu_argoneut_14_cs_lo = [ sqrt(x[5]**2.0+x[7]**2.0) for x in arr_raw ]

    # NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 3.6,1.5,1.5,0.28,0.01,0.01,0.03,0.03 ]
    ]

    arr_cs_cc_numubar_argoneut_14_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numubar_argoneut_14_energy_lo = [ x[1] for x in arr_raw ]
    arr_cs_cc_numubar_argoneut_14_energy_hi = [ x[2] for x in arr_raw ]
    arr_cs_cc_numubar_argoneut_14_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numubar_argoneut_14_cs_hi = [ sqrt(x[4]**2.0+x[6]**2.0) for x in arr_raw ]
    arr_cs_cc_numubar_argoneut_14_cs_lo = [ sqrt(x[5]**2.0+x[7]**2.0) for x in arr_raw ]
    #################################################################################


    #################################################################################
    # ArgoNeuT 2012
    # PRL 108, 161802 (2012)
    # https://arxiv.org/abs/1111.0103

    # NUMU
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 4.3,0.0,0.0,0.73,0.12,0.12 ]
    ]

    arr_cs_cc_numu_argoneut_12_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_argoneut_12_energy_lo = [ x[1] for x in arr_raw ]
    arr_cs_cc_numu_argoneut_12_energy_hi = [ x[2] for x in arr_raw ]
    arr_cs_cc_numu_argoneut_12_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_argoneut_12_cs_hi = [ x[4] for x in arr_raw ]
    arr_cs_cc_numu_argoneut_12_cs_lo = [ x[5] for x in arr_raw ]

    # No NUMUBAR data
    #################################################################################


    #################################################################################
    # GGM-SPS
    # PLB 104, 235 (1981)

    # NUMU
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 17.425742574257427,15.24752475247525,19.60396039603961,0.5795918367346939,0.6653061224489796,0.5020408163265309 ],
    [ 22.37623762376238,19.60396039603961,24.75247524752475,0.5632653061224491,0.6244897959183673,0.47755102040816344 ],
    [ 27.128712871287128,24.95049504950495,29.306930693069305,0.6367346938775511,0.7081632653061225,0.5489795918367351 ],
    [ 34.65346534653466,30.099009900990097,39.80198019801979,0.7142857142857144,0.8081632653061226,0.616326530612245 ],
    [ 44.75247524752476,40.00,49.60396039603961,0.6367346938775511,0.7346938775510203,0.536734693877551 ],
    [ 59.60396039603961,50.29702970297029,69.30693069306932,0.6489795918367347,0.7428571428571429,0.5469387755102044 ],
    [ 84.35643564356435,69.9009900990099,99.10891089108912,0.6326530612244898,0.7183673469387755,0.5326530612244897 ],
    [ 119.40594059405939,99.50495049504948,139.40594059405942,0.526530612244898,0.6122448979591837,0.4346938775510203 ]
    ]

    arr_cs_cc_numu_ggm_sps_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_ggm_sps_energy_lo = [ x[1] for x in arr_raw ]
    arr_cs_cc_numu_ggm_sps_energy_hi = [ x[2] for x in arr_raw ]
    arr_cs_cc_numu_ggm_sps_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_ggm_sps_cs_hi = [ x[4]-x[3] for x in arr_raw ]
    arr_cs_cc_numu_ggm_sps_cs_lo = [ x[3]-x[5] for x in arr_raw ]

    # NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 17.227722772277225,15.445544554455445,19.60396039603961,0.6257425742574321/2.,0.6891089108910933/2.,0.5465346534653568/2. ],
    [ 21.980198019801982,19.60396039603961,24.15841584158416,0.5702970297029792/2.,0.637623762376244/2.,0.49900990099011144/2. ],
    [ 27.524752475247524,24.95049504950495,29.306930693069305,0.5306930693069418/2.,0.5900990099009986/2.,0.4594059405940736/2. ],
    [ 34.45544554455445,29.900990099009903,39.20792079207921,0.6415841584158475/2.,0.7306930693069327/2.,0.5643564356435742/2. ],
    [ 44.75247524752474,40.00,49.801980198019805,0.5287128712871398/2.,0.6000000000000076/2.,0.4475247524752626/2. ],
    [ 59.999999999999986,50.0990099009901,69.50495049504951,0.6019801980198096/2.,0.6950495049504983/2.,0.5069306930693194/2. ],
    [ 84.85148514851485,69.9009900990099,99.60396039603958,0.629702970297036/2.,0.7366336633663373/2.,0.5069306930693189/2. ],
    [ 120.1980198019802,100.39603960396042,139.6039603960396,0.5702970297029792/2.0,0.6891089108910933/2.0,0.45742574257427204/2.0 ]
    ]

    arr_cs_cc_numubar_ggm_sps_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_sps_energy_lo = [ x[1] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_sps_energy_hi = [ x[2] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_sps_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_sps_cs_hi = [ x[4]-x[3] for x in arr_raw ]
    arr_cs_cc_numubar_ggm_sps_cs_lo = [ x[3]-x[5] for x in arr_raw ]
    #################################################################################


    #################################################################################
    # NuTeV
    # PRD 74, 012008 (2006)
    # https://arxiv.org/abs/hep-ex/0509010

    # NUMU
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 36.374169971923884,0.0,0.0,0.7030304380765633,0.0,0.0 ],
    [ 46.32002317393823,0.0,0.0,0.6851619947412988,0.0,0.0 ],
    [ 56.3093275101386,0.0,0.0,0.6816324256874193,0.0,0.0 ],
    [ 66.4534961451045,0.0,0.0,0.6792080752261689,0.0,0.0 ],
    [ 76.46396898257497,0.0,0.0,0.6826641115914258,0.0,0.0 ],
    [ 86.47221355675389,0.0,0.0,0.6853848210704578,0.0,0.0 ],
    [ 96.45483310307947,0.0,0.0,0.6796492713579038,0.0,0.0 ],
    [ 111.27724051873973,0.0,0.0,0.6710437185257812,0.0,0.0 ],
    [ 131.38286019876108,0.0,0.0,0.6558982129328401,0.6617808280226392,0.6496479343999287 ],
    [ 151.40380587370205,0.0,0.0,0.662810285663354,0.6690605641962655,0.656927670573555 ],
    [ 171.3311644903962,0.0,0.0,0.638838629172423,0.0,0.0 ],
    [ 191.57939302107934,0.0,0.0,0.6707540442978742,0.677737421453719,0.6648692009447836 ],
    [ 216.39333303623152,214.1205044788092,218.66616159365364,0.659354249298097,0.0,0.0 ],
    [ 246.48714292080754,244.06279245955702,248.91260751370353,0.6403115112081643,0.647294888364009,0.6336913409688488 ],
    [ 276.6154908864008,273.88809661749633,279.34177102366414,0.6826663398547169,0.6951691251838317,0.6716386648246355 ],
    [ 306.4251972013013,303.6966888007487,309.1514773385623,0.6198694237711126,0.6345759614956101,0.6047952226035026 ],
    [ 341.5927626008289,337.6509648380052,345.22706002941305,0.6251660056152233,0.6453874949864074,0.6038415259147019 ]
    ]

    arr_cs_cc_numu_nutev_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numu_nutev_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numu_nutev_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numu_nutev_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numu_nutev_cs_hi = [ x[4]-x[3] for x in arr_raw ]
    arr_cs_cc_numu_nutev_cs_lo = [ x[3]-x[5] for x in arr_raw ]

    # NUMUBAR
    # 'E [GEV] LOW','E [GEV] HIGH','SIG/E [10**-38CM**2/GEV]','error +','error -'

    arr_raw = [
    [ 35.156424083069666,0.0,0.0,0.3511742947546682,0.0,0.0 ],
    [ 45.279424216765456,0.0,0.0,0.34176433887428137,0.0,0.0 ],
    [ 55.42804937831454,0.0,0.0,0.3408106421854807,0.0,0.0 ],
    [ 65.29368510183163,0.0,0.0,0.3464704309461206,0.0,0.0 ],
    [ 75.43451134186014,0.0,0.0,0.34294309015553276,0.0,0.0 ],
    [ 85.4327287312269,0.0,0.0,0.34235482864655287,0.0,0.0 ],
    [ 95.58358215606755,0.0,0.0,0.3421364588439769,0.0,0.0 ],
    [ 110.43384286287268,0.0,0.0,0.3427224920896652,0.0,0.0 ],
    [ 130.26650028967427,0.0,0.0,0.337499442934177,0.0,0.0 ],
    [ 150.41445697223583,0.0,0.0,0.3363251481795087,0.0,0.0 ],
    [ 170.42871785730202,0.0,0.0,0.3410312402513481,0.0,0.0 ],
    [ 190.5555060385936,0.0,0.0,0.3328713400775434,0.0,0.0 ],
    [ 215.60229956771693,213.17572084317487,217.87401399349346,0.34831320468826577,0.0,0.0 ],
    [ 245.5423592851732,243.11912295556843,247.96782387806942,0.32853291144881647,0.3355185168679528,0.32191496947279263 ],
    [ 275.5247560051696,272.79513347297114,278.24992201078476,0.32272382904763997,0.3330206337180798,0.31243148090378337 ],
    [ 305.54169080618567,302.81429653727884,308.26908507509245,0.32831231338294903,0.34448950487989627,0.31139979499977677 ],
    [ 340.6747181246936,336.7351486251615,344.4627657203975,0.3222113284905741,0.3468447791791074,0.29831543295155727 ],
    ]

    arr_cs_cc_numubar_nutev_energy_central = [ x[0] for x in arr_raw ]
    arr_cs_cc_numubar_nutev_energy_lo = [ x[0]-x[1] for x in arr_raw ]
    arr_cs_cc_numubar_nutev_energy_hi = [ x[2]-x[0] for x in arr_raw ]
    arr_cs_cc_numubar_nutev_cs_central = [ x[3] for x in arr_raw ]
    arr_cs_cc_numubar_nutev_cs_hi = [ x[4]-x[3] for x in arr_raw ]
    arr_cs_cc_numubar_nutev_cs_lo = [ x[3]-x[5] for x in arr_raw ]
    #################################################################################



    #def Plot_Cross_Section_CC_HESE(output_format):
        
    #    print("Plot_Cross_Section_HESE: Plotting neutrino charged-current cross section data...")

    #    mpl.rcParams['xtick.labelsize']=26
    #    mpl.rcParams['ytick.labelsize']=26
    #    mpl.rcParams['legend.fontsize']=14
    #    mpl.rcParams['legend.borderpad']=0.8
    #    mpl.rcParams['ps.fonttype']=42
    #    mpl.rcParams['pdf.fonttype']=42
        # mpl.rcParams['legend.handlelength'] = 0
    #    mpl.rcParams['font.serif']='Palatino'

    #    fig = plt.figure(figsize=[9,9])
    #    ax = fig.add_subplot(1,1,1)

    #    ax.tick_params('both', length=10, width=2, which='major')
    #    ax.tick_params('both', length=5, width=1, which='minor')
    #    ax.tick_params(axis='both', which='major', pad=10)
    #    ax.tick_params(axis='y',which='minor', left='on')
    #    minorLocator = MultipleLocator(0.05)
    #    ax.yaxis.set_minor_locator(minorLocator)
     
        # IceCube HESE energy region
        # ax.fill_between( [2.5e4,5.e6], 
                         # [1.35,1.35],
                         # [0.0,0.0],
                         # facecolor=myLightBlue, alpha=0.95, edgecolor='None')

        # Horizontal dashed lines
    ax.plot( [ 0.1, 5.e6 ], [ 0.33, 0.33 ], '--', color='0.5' )
    ax.plot( [ 0.1, 5.e6 ], [ 0.675, 0.675 ], '--', color='0.5' )
        # ax.plot( [ 0.1, 5.e6 ], [ (0.33+0.675)/2.0, (0.33+0.675)/2.0 ], 'k--' )

        # DIS line (rough)
        # arr_dis_energy =[ 10.**x for x in np.linspace(-1.0,7.0,100) ]
        # avg_sigma_cc_low_energy = (0.33+0.675)/2.0
        # arr_dis_cs = [ avg_sigma_cc_low_energy*pow(energy/1.e4,0.5-1.0) for energy in arr_dis_energy ]
        # ax.plot( arr_dis_energy, arr_dis_cs, 'b-.', lw=1, dashes=(6,2))

        # DIS line (correct, but lacks the pdf contributions <x>; see Giunti & Kim)
        # avg_sigma_cc_low_energy = (0.33+0.675)/2.0
        # mass_nucleon = 938.91872965*1.e-3 # [GeV]
        # mass_w = 80.385 # [GeV] 
        # arr_dis_cs_raw = [ 2.*mass_nucleon * \
        #                 pow(1.0+(2.0*mass_nucleon*energy/mass_w/mass_w), -2.0) \
        #                 for energy in arr_dis_energy ] # sigma/E_nu
        # renorm_const = avg_sigma_cc_low_energy/arr_dis_cs_raw[0]
        # arr_dis_cs_renorm = [ x*renorm_const for x in arr_dis_cs_raw ]
        # ax.plot( arr_dis_energy, arr_dis_cs_renorm, 'b-.', lw=1, dashes=(6,2))

        # DIS line (full calculation)
    arr_dis_energy = [0.01, 0.012328467394420659, 0.015199110829529339, 0.018738174228603841, 0.023101297000831605, 0.028480358684358019, 0.035111917342151307, 0.043287612810830572, 0.0533669923120631, 0.065793322465756823, 0.081113083078968723, 0.10000000000000001, 0.12328467394420659, 0.15199110829529339, 0.18738174228603841, 0.23101297000831603, 0.2848035868435802, 0.35111917342151311, 0.43287612810830595, 0.533669923120631, 0.65793322465756821, 0.81113083078968728, 1.0, 1.2328467394420659, 1.5199110829529332, 1.873817422860385, 2.3101297000831602, 2.8480358684358018, 3.5111917342151311, 4.3287612810830574, 5.3366992312063131, 6.5793322465756825, 8.1113083078968717, 10.0, 12.32846739442066, 15.199110829529348, 18.73817422860385, 23.101297000831604, 28.48035868435802, 35.111917342151308, 43.287612810830616, 53.366992312063125, 65.79332246575683, 81.113083078968728, 100.0, 123.28467394420659, 151.99110829529332, 187.3817422860383, 231.01297000831579, 284.80358684358049, 351.11917342151344, 432.87612810830615, 533.66992312063121, 657.93322465756819, 811.13083078968725, 1000.0, 1232.8467394420659, 1519.9110829529332, 1873.817422860383, 2310.1297000831628, 2848.0358684358048, 3511.1917342151346, 4328.7612810830615, 5336.699231206313, 6579.3322465756828, 8111.3083078968721, 10000.0, 12328.467394420659, 15199.110829529331, 18738.174228603832, 23101.297000831626, 28480.358684358049, 35111.917342151348, 43287.612810830615, 53366.992312063128, 65793.322465756821, 81113.083078968717, 100000.0, 123284.67394420659, 151991.10829529332, 187381.74228603867, 231012.97000831628, 284803.58684358047, 351119.17342151346, 432876.12810830615, 533669.92312063125, 657933.22465756827, 811130.83078968723, 1000000.0, 1232846.7394420684, 1519911.0829529332, 1873817.4228603868, 2310129.700083158, 2848035.8684358047, 3511191.7342151273, 4328761.2810830614, 5336699.2312063016, 6579332.2465756824, 8111308.3078968888, 10000000.0] # [GeV]
    arr_dis_cs_raw = [3.9146146595649999e-37, 3.1906829797467376e-37, 2.6034796683606003e-37, 2.1271809505327561e-37, 1.7408403659533781e-37, 1.4274675944282592e-37, 1.1732812628873826e-37, 9.6710287408230948e-38, 7.9986520343624425e-38, 6.6421354454231947e-38, 5.5418226704119274e-38, 4.6493246269650335e-38, 3.9253914196547486e-38, 3.3381862249764967e-38, 2.8618676171441036e-38, 2.4751494585265832e-38, 2.1596794867858922e-38, 1.8999376590565093e-38, 1.6841797050259346e-38, 1.5037139966395841e-38, 1.352049006634688e-38, 1.2242568022959765e-38, 1.1166638957990549e-38, 1.0260347696355845e-38, 9.4947788121194603e-39, 8.8482957604583784e-39, 8.3019307547061459e-39, 7.8394611601493065e-39, 7.4477312241119139e-39, 7.1147416417711783e-39, 6.829776416483251e-39, 6.5899901542038582e-39, 6.3863276747282906e-39, 6.2103953549528084e-39, 6.0597141824432318e-39, 5.9308560921854949e-39, 5.8196213023870665e-39, 5.7243577413567575e-39, 5.6396491467004515e-39, 5.5655877780223423e-39, 5.5013823742379428e-39, 5.4447810649852793e-39, 5.3941244978134042e-39, 5.3480094251799645e-39, 5.3058749244755088e-39, 5.2656042172763661e-39, 5.2295640164902756e-39, 5.1948058874061e-39, 5.161662376476451e-39, 5.1292258435352807e-39, 5.0963766155801884e-39, 5.0632762725742924e-39, 5.0285774462847825e-39, 4.990827944314087e-39, 4.9502074480457819e-39, 4.9051839504853769e-39, 4.8548288265854508e-39, 4.7987460651454269e-39, 4.7356405245702886e-39, 4.6640006715605618e-39, 4.5828965669436114e-39, 4.491392600588347e-39, 4.3886793018736097e-39, 4.2754258372420212e-39, 4.1496470206296528e-39, 4.0127270744319147e-39, 3.865228641804291e-39, 3.7081277878136196e-39, 3.5428258076465466e-39, 3.3706297729674556e-39, 3.1927363445945573e-39, 3.0112344464795341e-39, 2.8281442010715574e-39, 2.6462902599409864e-39, 2.4674392089724955e-39, 2.2923996637739577e-39, 2.1224337035251038e-39, 1.9585813760504172e-39, 1.8018021614928074e-39, 1.6526534084105611e-39, 1.5116647989100753e-39, 1.379065946196128e-39, 1.254925926257109e-39, 1.139220521117625e-39, 1.0318351952295794e-39, 9.3250009423533595e-40, 8.4086229607135346e-40, 7.5656016531790672e-40, 6.7917639646467664e-40, 6.0823935419015025e-40, 5.4322249577114891e-40, 4.8361149226985023e-40, 4.2902302294424997e-40, 3.7907610791368109e-40, 3.3343938131109214e-40, 2.9186004033953091e-40, 2.5417013450426419e-40, 2.2017812039912533e-40, 1.8966629934023095e-40, 1.6246463567699999e-40] # [cm^2 GeV^-1]
    avg_sigma_cc_low_energy = (0.33+0.675)/2.0
        # renorm_const = avg_sigma_cc_low_energy/(arr_dis_cs_raw[44]*1.e38)
    renorm_const = 1.0
        # print arr_dis_energy[44]
    arr_dis_cs_renorm = [ x*renorm_const*1.e38 for x in arr_dis_cs_raw ]
        # print arr_dis_cs_renorm
    ax.plot( arr_dis_energy[34:len(arr_dis_energy)], arr_dis_cs_renorm[34:len(arr_dis_cs_renorm)], 'k-.', lw=1, dashes=(6,2))

        
        # NUMU

        # ax.plot( arr_cs_cc_numu_t2k_c_2015_energy_central, arr_cs_cc_numu_t2k_c_2015_cs_central, linestyle='None', marker='+', color='0.5', mew=4, markersize=11, label='T2K (C)' )         
    
    col_lines = []
    cosmic_lines = []
        
    ax.errorbar( arr_cs_cc_numu_t2k_fe_2014_energy_central, arr_cs_cc_numu_t2k_fe_2014_cs_central, xerr=[arr_cs_cc_numu_t2k_fe_2014_energy_lo, arr_cs_cc_numu_t2k_fe_2014_energy_hi], yerr=[arr_cs_cc_numu_t2k_fe_2014_cs_lo, arr_cs_cc_numu_t2k_fe_2014_cs_hi], linestyle='None', marker='+', color='k', mew=1, markersize=13 ) # , label='T2K (Fe)'
    ax.plot( arr_cs_cc_numu_t2k_fe_2014_energy_central, arr_cs_cc_numu_t2k_fe_2014_cs_central, linestyle='None', marker='+', color='k', mew=3, markersize=13 ) # , label='T2K (Fe)'
    col_line = ax.plot( [1.e-2, -1.0], linestyle='None', marker='+', color='k', mew=3, markersize=13, label='T2K (Fe) 14')
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_t2k_ch_2014_energy_central, arr_cs_cc_numu_t2k_ch_2014_cs_central, xerr=[arr_cs_cc_numu_t2k_ch_2014_energy_lo, arr_cs_cc_numu_t2k_ch_2014_energy_hi], yerr=[arr_cs_cc_numu_t2k_ch_2014_cs_lo, arr_cs_cc_numu_t2k_ch_2014_cs_hi], linestyle='None', marker='+', color='0.5', mew=1, markersize=13 ) # , label='T2K (CH)'
    ax.plot( arr_cs_cc_numu_t2k_ch_2014_energy_central, arr_cs_cc_numu_t2k_ch_2014_cs_central, linestyle='None', marker='+', color='0.5', mew=3, markersize=13 ) # , label='T2K (Fe)'
    col_line = ax.plot( [1.e-2, -1.0], linestyle='None', marker='+', color='0.5', mew=3, markersize=13, label='T2K (CH) 14')
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_t2k_c_2013_energy_central, arr_cs_cc_numu_t2k_c_2013_cs_central, xerr=[arr_cs_cc_numu_t2k_c_2013_energy_lo, arr_cs_cc_numu_t2k_c_2013_energy_hi], yerr=[arr_cs_cc_numu_t2k_c_2013_cs_lo, arr_cs_cc_numu_t2k_c_2013_cs_hi], linestyle='None', marker='*', color='k', mew=1, markerfacecolor='w', markersize=13 ) # , label='T2K (C)'
    col_line = ax.plot( [1.e-2, -1.0], linestyle='None', marker='*', color='k', mew=1, markerfacecolor='w', markersize=13, label='T2K (C) 13')
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_argoneut_14_energy_central, arr_cs_cc_numu_argoneut_14_cs_central, xerr=[arr_cs_cc_numu_argoneut_14_energy_lo, arr_cs_cc_numu_argoneut_14_energy_hi], yerr=[arr_cs_cc_numu_argoneut_14_cs_lo, arr_cs_cc_numu_argoneut_14_cs_hi], linestyle='None', color='k', mew=1, markerfacecolor='w', markersize=10, fmt='^' ) # , label='ArgoNeuT 2014'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', mew=1, markerfacecolor='w', markersize=10, marker=r'^', label='ArgoNeuT 14' )
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_argoneut_12_energy_central, arr_cs_cc_numu_argoneut_12_cs_central, yerr=[arr_cs_cc_numu_argoneut_12_cs_lo, arr_cs_cc_numu_argoneut_12_cs_hi], linestyle='None', color='0.5', markeredgecolor='0.5',markersize=8, fmt='o', zorder=10 ) # , label='ArgoNeuT 2012'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='0.5', markeredgecolor='0.5', markersize=8, marker=r'o', label='ArgoNeuT 12' )    
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_anl_energy_central, arr_cs_cc_numu_anl_cs_central, yerr=[arr_cs_cc_numu_anl_cs_lo, arr_cs_cc_numu_anl_cs_hi], linestyle='None', color='k', markersize=13, fmt='*' ) # , label='ANL'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', marker='*', color='k', markersize=13, label='ANL 79' )
    col_lines.append(col_line)
            
        # ax.errorbar( arr_cs_cc_numu_bebc_1_energy_central, arr_cs_cc_numu_bebc_1_cs_central, yerr=[arr_cs_cc_numu_bebc_1_cs_lo, arr_cs_cc_numu_bebc_1_cs_hi], linestyle='None', color='k', markersize=13, label='BEBC', fmt='--', marker=r'$\bigcirc$' ) 

    ax.errorbar( arr_cs_cc_numu_bebc_2_energy_central, arr_cs_cc_numu_bebc_2_cs_central, yerr=[arr_cs_cc_numu_bebc_2_cs_lo, arr_cs_cc_numu_bebc_2_cs_hi], linestyle='None', color='k', markersize=10, fmt='', marker=r'$\bigcirc$' ) # , label='BEBC'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', markersize=10, marker=r'$\bigcirc$', label='BEBC 79' )
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_bnl_energy_central, arr_cs_cc_numu_bnl_cs_central, yerr=[arr_cs_cc_numu_bnl_cs_lo, arr_cs_cc_numu_bnl_cs_hi], linestyle='None', color='k', markeredgecolor='k', markersize=10, fmt='^' ) # , label='BNL'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', markeredgecolor='k', markersize=10, marker=r'^', label='BNL 82' )    
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_ccfr_energy_central, arr_cs_cc_numu_ccfr_cs_central, yerr=[arr_cs_cc_numu_ccfr_cs_lo, arr_cs_cc_numu_ccfr_cs_hi], linestyle='None', color='k', mew=1, markerfacecolor='w', markersize=8, fmt='', marker='d', zorder= 2 ) # , label='CCFR'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', mew=1, markerfacecolor='w', markersize=8, marker='d' , label='CCFR 97' )
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_cdhs_energy_central, arr_cs_cc_numu_cdhs_cs_central, yerr=[arr_cs_cc_numu_cdhs_cs_lo, arr_cs_cc_numu_cdhs_cs_hi], linestyle='None', color='k', mew=1, markerfacecolor='w', markersize=8, fmt='s' ) # , label='CDHS'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', mew=1, markerfacecolor='w', markersize=8, marker=r's', label='CDHS 87' )
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_ggm_sps_energy_central, arr_cs_cc_numu_ggm_sps_cs_central, yerr=[arr_cs_cc_numu_ggm_sps_cs_lo, arr_cs_cc_numu_ggm_sps_cs_hi], linestyle='None', color='0.5', markeredgecolor='0.5', markersize=8, fmt='--s' ) # , label='GGM-SPS'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='0.5', markeredgecolor='0.5', markersize=8, marker=r's', label='GGM-SPS 81' )
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_ggm_ps_energy_central, arr_cs_cc_numu_ggm_ps_cs_central, yerr=[arr_cs_cc_numu_ggm_ps_cs_lo, arr_cs_cc_numu_ggm_ps_cs_hi], linestyle='None', color='k', markersize=8, fmt='s' ) # , label='GGM-PS'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', markersize=8, marker=r's', label='GGM-PS 79' )
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_ihep_itep_energy_central, arr_cs_cc_numu_ihep_itep_cs_central, yerr=[arr_cs_cc_numu_ihep_itep_cs_lo, arr_cs_cc_numu_ihep_itep_cs_hi], linestyle='None', color='k', markersize=10, fmt='v' ) # , label='IHEP-ITEP'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', markersize=10, marker=r'v', label='IHEP-ITEP 79' )
    col_lines.append(col_line)
            
    ax.errorbar( arr_cs_cc_numu_ihep_jinr_energy_central, arr_cs_cc_numu_ihep_jinr_cs_central, yerr=[arr_cs_cc_numu_ihep_jinr_cs_lo, arr_cs_cc_numu_ihep_jinr_cs_hi], linestyle='None', color='0.5', markeredgecolor='0.5', markersize=10, fmt='v' ) # , label='IHEP-JINR'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='0.5', markeredgecolor='0.5', markersize=10, marker=r'v', label='IHEP-JINR 96' )
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_minos_energy_central, arr_cs_cc_numu_minos_cs_central, yerr=[arr_cs_cc_numu_minos_cs_lo, arr_cs_cc_numu_minos_cs_hi], linestyle='None', color='k', markersize=8, fmt='o', zorder=10 ) # , label='MINOS'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', markersize=8, marker=r'o', label='MINOS 10' )
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_nomad_energy_central, arr_cs_cc_numu_nomad_cs_central, xerr=[arr_cs_cc_numu_nomad_energy_lo, arr_cs_cc_numu_nomad_energy_hi], yerr=[arr_cs_cc_numu_nomad_cs_lo, arr_cs_cc_numu_nomad_cs_hi], linestyle='None', color='0.5', markeredgecolor='0.5', markersize=10, fmt='^' ) # , label='NOMAD'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='0.5', markeredgecolor='0.5', markersize=10, marker=r'^', label='NOMAD 08' )    
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_nutev_energy_central, arr_cs_cc_numu_nutev_cs_central, xerr=[arr_cs_cc_numu_nutev_energy_lo, arr_cs_cc_numu_nutev_energy_hi], yerr=[arr_cs_cc_numu_nutev_cs_lo, arr_cs_cc_numu_nutev_cs_hi], linestyle='None', color='k', mew=1, markerfacecolor='k', markersize=8, fmt='', marker='d', zorder=0 ) # , label='NuTeV'
    col_line = ax.plot( [1.e-2, -1.0], linestyle='None', color='k', mew=1, markerfacecolor='k', markersize=8, marker='d', label='NuTeV 06')
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_sciboone_energy_central, arr_cs_cc_numu_sciboone_cs_central, xerr=[arr_cs_cc_numu_sciboone_energy_lo, arr_cs_cc_numu_sciboone_energy_hi], yerr=[arr_cs_cc_numu_sciboone_cs_lo, arr_cs_cc_numu_sciboone_cs_hi], linestyle='None', color='k', markersize=10, fmt='x' ) # , label='SciBooNE'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', markersize=10, marker=r'x', label='SciBooNE 11' )    
    col_lines.append(col_line)
    
    ax.errorbar( arr_cs_cc_numu_skat_energy_central, arr_cs_cc_numu_skat_cs_central, yerr=[arr_cs_cc_numu_skat_cs_lo, arr_cs_cc_numu_skat_cs_hi], linestyle='None', color='k', markersize=11, fmt='', marker=r'$\otimes$' ) # , label='SKAT'
    col_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', markersize=11, marker=r'$\otimes$', label='SKAT 79' )    
    col_lines.append(col_line)
    
        # ax.errorbar( arr_cs_cc_numu_ic_energy_central, arr_cs_cc_numu_ic_cs_central, xerr=[arr_cs_cc_numu_ic_energy_lo, arr_cs_cc_numu_ic_energy_hi], yerr=[arr_cs_cc_numu_ic_cs_lo, arr_cs_cc_numu_ic_cs_hi], linestyle='None', color='b', markeredgecolor='b', markersize=10, fmt='--', marker=r'$\ast$' ) # , label='IceCube'
        # ax.plot( [ 1.e-2, -1.0], linestyle='None', color='b', markeredgecolor='b', markersize=10, marker=r'$\ast$', label='IC HESE sh. 17' )

        # ax.errorbar( arr_cs_cc_nuavg_ic_energy_central[1:4], arr_cs_cc_nuavg_ic_cs_central[1:4], xerr=[arr_cs_cc_nuavg_ic_energy_lo[1:4], arr_cs_cc_nuavg_ic_energy_hi[1:4]], yerr=[arr_cs_cc_nuavg_ic_cs_lo[1:4], arr_cs_cc_nuavg_ic_cs_hi[1:4]], linestyle='None', color='b', markeredgecolor='b', markersize=10, fmt='--', marker=r'$\ast$' ) # , label='IceCube'
        # ax.errorbar( [arr_cs_cc_nuavg_ic_energy_central[4]], [arr_cs_cc_nuavg_ic_cs_central[4]-arr_cs_cc_nuavg_ic_cs_lo[4]], xerr=[[arr_cs_cc_nuavg_ic_energy_lo[4]], [arr_cs_cc_nuavg_ic_energy_hi[3]]], linestyle='None', color='b', markeredgecolor='b', markersize=10, fmt='--', marker=None ) # , label='IceCube'
        # ax.plot( [ 1.e-2, -1.0], linestyle='None', color='b', markeredgecolor='b', markersize=10, marker=r'$\ast$', label=r'IC HESE showers 17 (avg. of $\nu$, $\bar{\nu}$)' )
        # ax.plot( [arr_cs_cc_nuavg_ic_energy_central[4], arr_cs_cc_nuavg_ic_energy_central[4]], [arr_cs_cc_nuavg_ic_cs_central[4]-arr_cs_cc_nuavg_ic_cs_lo[4], 0.2*arr_cs_cc_nuavg_ic_cs_central[4]], 'b-' )
        # ax.annotate( r'$>$', xy = (1.029*arr_cs_cc_nuavg_ic_energy_central[4], 0.18*arr_cs_cc_nuavg_ic_cs_central[4]), xycoords='data', color='b', fontsize=20, horizontalalignment='center', rotation=90 )
        # Uncomment the next two lines if using lo-res points
        # ax.plot( [arr_cs_cc_nuavg_ic_energy_central[3], arr_cs_cc_nuavg_ic_energy_central[3]], [arr_cs_cc_nuavg_ic_cs_central[3]-arr_cs_cc_nuavg_ic_cs_lo[3], 0.8*arr_cs_cc_nuavg_ic_cs_central[3]], 'b-' )
        # ax.annotate( r'$>$', xy = (1.025*arr_cs_cc_nuavg_ic_energy_central[3], 0.721*arr_cs_cc_nuavg_ic_cs_central[3]), xycoords='data', color='b', fontsize=20, horizontalalignment='center', rotation=90 )

    ax.errorbar( arr_cs_cc_nuavg_ic_energy_central[0:3], arr_cs_cc_nuavg_ic_cs_central[0:3], xerr=[arr_cs_cc_nuavg_ic_energy_lo[0:3], arr_cs_cc_nuavg_ic_energy_hi[0:3]], yerr=[arr_cs_cc_nuavg_ic_cs_lo[0:3], arr_cs_cc_nuavg_ic_cs_hi[0:3]], linestyle='None', color='k', markeredgecolor='k', markersize=10, fmt='', marker=r'$\ast$' ) # , label='IceCube'
    ax.errorbar( [arr_cs_cc_nuavg_ic_energy_central[3]], [arr_cs_cc_nuavg_ic_cs_central[3]-arr_cs_cc_nuavg_ic_cs_lo[3]], xerr=[[arr_cs_cc_nuavg_ic_energy_lo[3]], [arr_cs_cc_nuavg_ic_energy_hi[3]]], linestyle='None', color='k', markeredgecolor='k', markersize=10, fmt='', marker=None ) # , label='IceCube'
    ax.plot( [arr_cs_cc_nuavg_ic_energy_central[3], arr_cs_cc_nuavg_ic_energy_central[3]], [arr_cs_cc_nuavg_ic_cs_central[3]-arr_cs_cc_nuavg_ic_cs_lo[3], 0.168*arr_cs_cc_nuavg_ic_cs_central[3]], 'k-' )
    ax.annotate( r'$>$', xy = (1.025*arr_cs_cc_nuavg_ic_energy_central[3], 0.151*arr_cs_cc_nuavg_ic_cs_central[3]), xycoords='data', color='k', fontsize=20, horizontalalignment='center', rotation=90 )
    cosmic_line = ax.plot( [ 1.e-2, -1.0], linestyle='None', color='k', markeredgecolor='k', markersize=10, marker=r'$\ast$', label=r'IC HESE showers 17 (avg. of $\nu$, $\bar{\nu}$)' )
    
    cosmic_lines.append(cosmic_line)


        # NUMUBAR

        # ax.errorbar( arr_cs_cc_numubar_bebc_1_energy_central, arr_cs_cc_numubar_bebc_1_cs_central, yerr=[arr_cs_cc_numubar_bebc_1_cs_lo, arr_cs_cc_numubar_bebc_1_cs_hi], linestyle='None', color='k', markersize=13, fmt='--', marker=r'$\bigcirc$' ) # label='BEBC', 

    ax.errorbar( arr_cs_cc_numubar_argoneut_14_energy_central, arr_cs_cc_numubar_argoneut_14_cs_central, xerr=[arr_cs_cc_numubar_argoneut_14_energy_lo, arr_cs_cc_numubar_argoneut_14_energy_hi], yerr=[arr_cs_cc_numubar_argoneut_14_cs_lo, arr_cs_cc_numubar_argoneut_14_cs_hi], linestyle='None', color='k', mew=1, markerfacecolor='w', markersize=10, fmt='^' ) # , label='ArgoNeuT 2014'

    ax.errorbar( arr_cs_cc_numubar_ccfr_energy_central, arr_cs_cc_numubar_ccfr_cs_central, yerr=[arr_cs_cc_numubar_ccfr_cs_lo, arr_cs_cc_numubar_ccfr_cs_hi], linestyle='None', color='k', mew=1, markerfacecolor='w', markersize=8, fmt='', marker='d', zorder=2 ) # , label='CCFR'

    ax.errorbar( arr_cs_cc_numubar_cdhs_energy_central, arr_cs_cc_numubar_cdhs_cs_central, yerr=[arr_cs_cc_numubar_cdhs_cs_lo, arr_cs_cc_numubar_cdhs_cs_hi], linestyle='None', color='k', mew=1, markerfacecolor='w', markersize=8, fmt='--s' ) # , label='CDHS'

    ax.errorbar( arr_cs_cc_numubar_ggm_sps_energy_central, arr_cs_cc_numubar_ggm_sps_cs_central, yerr=[arr_cs_cc_numubar_ggm_sps_cs_lo, arr_cs_cc_numubar_ggm_sps_cs_hi], linestyle='None', color='0.5', markeredgecolor='0.5', markersize=8, fmt='s' ) # , label='GGM-SPS'

    ax.errorbar( arr_cs_cc_numubar_ggm_ps_energy_central, arr_cs_cc_numubar_ggm_ps_cs_central, yerr=[arr_cs_cc_numubar_ggm_ps_cs_lo, arr_cs_cc_numubar_ggm_ps_cs_hi], linestyle='None', color='k', markersize=8, fmt='s' ) # label='GGM-PS',
        
    ax.errorbar( arr_cs_cc_numubar_ihep_itep_energy_central, arr_cs_cc_numubar_ihep_itep_cs_central, yerr=[arr_cs_cc_numu_ihep_itep_cs_lo, arr_cs_cc_numu_ihep_itep_cs_hi], linestyle='None', color='k', markersize=10, fmt='v' ) # label='IHEP-ITEP',

    ax.errorbar( arr_cs_cc_numubar_ihep_jinr_energy_central, arr_cs_cc_numubar_ihep_jinr_cs_central, yerr=[arr_cs_cc_numubar_ihep_jinr_cs_lo, arr_cs_cc_numubar_ihep_jinr_cs_hi], linestyle='None', color='0.5', markeredgecolor='0.5', markersize=10, fmt='v' ) # label='IHEP-JINR',

    ax.errorbar( arr_cs_cc_numubar_minos_energy_central, arr_cs_cc_numubar_minos_cs_central, yerr=[arr_cs_cc_numubar_minos_cs_lo, arr_cs_cc_numubar_minos_cs_hi], linestyle='None', color='k', markersize=8, fmt='o' ) # label='MINOS',

    ax.errorbar( arr_cs_cc_numubar_nutev_energy_central, arr_cs_cc_numubar_nutev_cs_central, xerr=[arr_cs_cc_numubar_nutev_energy_lo, arr_cs_cc_numubar_nutev_energy_hi], yerr=[arr_cs_cc_numubar_nutev_cs_lo, arr_cs_cc_numubar_nutev_cs_hi], linestyle='None', color='k', mew=1, markerfacecolor='k', markersize=8, fmt='--', marker='d', zorder=0 ) # , label='NuTeV'

        # ax.errorbar( arr_cs_cc_numubar_ic_energy_central, arr_cs_cc_numubar_ic_cs_central, xerr=[arr_cs_cc_numubar_ic_energy_lo, arr_cs_cc_numubar_ic_energy_hi], yerr=[arr_cs_cc_numubar_ic_cs_lo, arr_cs_cc_numubar_ic_cs_hi], linestyle='None', color='b', markeredgecolor='b', markersize=10, fmt='--', marker=r'$\ast$' ) # label='IceCube', 

        # ax.legend(loc='upper right', numpoints=1, ncol=2, columnspacing=0.5, frameon=False) 

        # Uncomment the next line for log y-axis
    #    leg = ax.legend(loc='lower left', numpoints=1, ncol=2, columnspacing=0.5, frameon=False)
    #    setp(leg.get_texts()[19], color='b')

        # # get handles
        # handles, labels = ax.get_legend_handles_labels()
        # # remove the errorbars
        # handles = [h[0] for h in handles]
        # # use them in the legend
        # ax.legend(handles, labels, loc='upper right',numpoints=1)

        # handles, labels = ax.get_legend_handles_labels()
        # for h in handles: h.set_linestyle("")
        # ax.legend(handles, labels)

    #    pylab.xlim([0.1, 5.e6])
    #    ax.set_xscale('log')

        # pylab.ylim([0.0, 1.35])

        # Uncomment the next two lines for log y-axis
    #    pylab.ylim([0.01, 2.0])
    #    ax.set_yscale('log')

        # xlabels = [item.get_text() for item in ax.get_xticklabels()]
    #    xlabels = [r'$10^{-2}$',r'$10^{-1}$',1,r'$10^1$',r'$10^2$',r'$10^3$',r'$10^4$',r'$10^5$',r'$10^6$']
    #    ax.set_xticklabels(xlabels)
     
    #    ax.set_xlabel(r'Neutrino energy $E_\nu$ [GeV]', fontsize=25)
    #    ax.set_ylabel('r$\sigma_{\nu N}^{\rm CC} / E_\nu$ [$10^{-38}$ cm$^2$\, GeV$^{-1}$\, nucleon$^{-1}$]', fontsize=25)
     
        # ax.annotate( r'$\nu$', xy = (0.30,0.75), xycoords = 'axes fraction', color='k', fontsize=24, horizontalalignment='left' )
        # ax.annotate( r'$\bar{\nu}$', xy = (0.30,0.10), xycoords = 'axes fraction', color='k', fontsize=24, horizontalalignment='left' )

        # Uncomment the next two lines for log y-axis
#     ax.annotate( r'$\nu$', xy = (0.30,0.89), xycoords = 'axes fraction', color='k', fontsize=24, horizontalalignment='left' )
#     ax.annotate( r'$\bar{\nu}$', xy = (0.30,0.52), xycoords = 'axes fraction', color='k', fontsize=24, horizontalalignment='left' )

        # ax.annotate( r'Avg. of $\nu$, $\bar{\nu}$ DIS', xy = (0.45,0.42), xycoords = 'axes fraction', color='b', fontsize=20, horizontalalignment='left', rotation = -30. )
        # ax.annotate( r'Avg. of $\nu$, $\bar{\nu}$ DIS', xy = (0.45,0.42), xycoords = 'axes fraction', color='b', fontsize=20, horizontalalignment='left')
        # ax.annotate( r'$\propto E_\nu^{-0.5}$', xy = (0.685,0.29), xycoords = 'axes fraction', color='b', fontsize=20, horizontalalignment='left' )
        # ax.annotate( r'Avg. of $\nu$, $\bar{\nu}$', xy = (0.43,0.10), xycoords = 'axes fraction', color='b', fontsize=24, horizontalalignment='left')

#     ax.annotate( r'DIS', xy = (0.57,0.73), xycoords = 'axes fraction', color='k', fontsize=18, horizontalalignment='left', rotation=-15)

        # pylab.savefig(os.getcwd()+"/"+"cross_sections_6yr."+output_format,bbox_inches='tight')
        # pylab.savefig(os.getcwd()+"/"+"cross_sections_talk."+output_format,bbox_inches='tight')
        # plt.close()

    #    return 0
    return col_lines, cosmic_lines


    #Plot_Cross_Section_CC_HESE('pdf')

