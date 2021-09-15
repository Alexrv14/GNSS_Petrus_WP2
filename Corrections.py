#!/usr/bin/env python

########################################################################
# PETRUS/SRC/Corrections.py:
# This is the Corrections Module of PETRUS tool
#
# Project:        PETRUS
# File:           Corrections.py
# Date(YY/MM/DD): 16/02/21
#
# Author: GNSS Academy
# Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################

# Import External and Internal functions and Libraries
#----------------------------------------------------------------------
import sys, os
# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)
from collections import OrderedDict
from COMMON import GnssConstants as Const
from COMMON import Iono, Tropo
from InputOutput import RcvrIdx, SatIdx, LosIdx
import numpy as np

def runCorrectMeas(Conf, Rcvr, PreproObsInfo, SatInfo, LosInfo):

    # Purpose: correct GNSS preprocessed measurements and compute
    #          pseudo range residuals

    #          More in detail, this function handles the following:
    #          tasks:

    #             *  Correct the satellite navigation position and clock using EGNOS Fast-Long-Term (FLT) corrections: FC and LTC.
    #             *  Estimate the Slant Ionospheric Delay (UISD) using MOPS guidelines interpolation criteria for IGP Selection
    #             *  Estimate the Slant Troposphere delay (STD) using MOPS model (ZTD) and its mapping function. 
    #             *  Correct the Pre-processed measurements from Geometrical Range, Satellite clock, ionosphere and troposphere. 
    #             *  Build the Corrected Measurements and Measurement Residuals
    #             *  Estimate all Range level Sigmas to build the Sigma UERE:
    #                   -  Estimate the SigmaUIRE from MT26 information
    #                   -  Estimate the SigmaFLT from UDRE and MT28 
    #                   -  Estimate the SigmaTRO budget in line with MOPS.
    #                   -  Estimate the SigmaAirborne budget in line with MOPS 
    #             *  Estimate the Sigma UERE budget in line with MOPS

    # Parameters
    # ==========
    # Conf: dict
    #       Configuration dictionary
    # Rcvr: list
    #       Receiver information: position, masking angle...
    # PreproObsInfo: dict
    #                Preprocessed observations for current epoch per sat
    #                PreproObsInfo["G01"]["C1"]
    # SatInfo: dict
    #          dictionary containing the split lines of the SAT file
    #          SatInfo["G01"][1] is the second field of the line
    #          containing G01 info
    # LosInfo: dict
    #          dictionary containing the split lines of the LOS file
    #          SatInfo["G01"][1] is the second field of the line
    #          containing G01 info

    # Returns
    # =======
    # CorrInfo: dict
    #           Corrected measurements for current epoch per sat
    #           CorrInfo["G01"]["CorrectedPsr"]

    # Initialize output
    CorrInfo = OrderedDict({})

    # Initialize some values
    ResSum = 0.0
    ResN = 0
    EntGpsSum = 0.0
    EntGpsN = 0

    # Loop over satellites
    for SatLabel, SatPrepro in PreproObsInfo.items():
        # If satellite is in convergence
        if(SatPrepro["Status"] == 1):
            # Initialize output info
            SatCorrInfo = {
                "Sod": 0.0,             # Second of day
                "Doy": 0,               # Day of year
                "Elevation": 0.0,       # Elevation
                "Azimuth": 0.0,         # Azimuth
                "IppLon": 0.0,          # IPP Longitude
                "IppLat": 0.0,          # IPP Latitude
                "Flag": 1,              # 0: Not Used 1: Used for PA 2: Used for NPA
                "SatX": 0.0,            # X-Component of the Satellite Position corrected with SBAS LTC
                "SatY": 0.0,            # Y-Component of the Satellite Position corrected with SBAS LTC
                "SatZ": 0.0,            # Z-Component of the Satellite Position corrected with SBAS LTC
                "SatClk": 0.0,          # Satellite Clock corrected with SBAS FLT
                "Uisd": 0.0,            # User Ionospheric Slant Delay
                "Std": 0.0,             # Slant Tropospheric Delay
                "CorrPsr": 0.0,         # Pseudo Range corrected from delays
                "GeomRange": 0.0,       # Geometrical Range (distance between Satellite Position and Receiver Reference Position)
                "PsrResidual": 0.0,     # Pseudo Range Residual
                "RcvrClk": 0.0,         # Receiver Clock estimation
                "SigmaFlt": 0,          # Sigma of the residual error associated to the fast and long-term correction (FLT)
                "SigmaUire": 0,         # User Ionospheric Range Error Sigma
                "SigmaTropo": 0,        # Sigma of the Tropospheric error 
                "SigmaAirborne": 0.0,   # Sigma Airborne Error
                "SigmaNoiseDiv": 0.0,   # Sigma of the receiver noise + divergence
                "SigmaMultipath": 0.0,  # Sigma of the receiver multipath
                "SigmaUere": 0.0,       # Sigma User Equivalent Range Error (Sigma of the total residual error associated to the satellite)
                "EntGps": 0.0,          # ENT to GPS Offset

            } # End of SatCorrInfo

            # Prepare outputs
            # Get SoD
            SatCorrInfo["Sod"] = SatPrepro["Sod"]
            # Get DoY
            SatCorrInfo["Doy"] = SatPrepro["Doy"]
            # Get Elevation
            SatCorrInfo["Elevation"] = SatPrepro["Elevation"]
            # Get Azimuth
            SatCorrInfo["Azimuth"] = SatPrepro["Azimuth"]

            # If SBAS information is available for current satellite
            if (SatLabel in SatInfo) and (SatLabel in LosInfo):
                # Get IPP Longitude
                SatCorrInfo["IppLon"] = float(LosInfo[SatLabel][LosIdx["IPPLON"]])
                # Get IPP Latitude
                SatCorrInfo["IppLat"] = float(LosInfo[SatLabel][LosIdx["IPPLAT"]])
                # Get Monitoring Status from UDREi
                if int(SatInfo[SatLabel][SatIdx["UDREI"]]) < 14:             
                    if int(SatInfo[SatLabel][SatIdx["UDREI"]]) < 12:
                        SatCorrInfo["Flag"] = 1                         # Satellite Monitored and Usable for PA
                    else:
                        SatCorrInfo["Flag"] = 2                         # Satellite Monitored and Usable for NPA
                elif int(SatInfo[SatLabel][SatIdx["UDREI"]]) == 14:          
                    SatCorrInfo["Flag"] = 0                             # Satellite Non-Monitored
                elif int(SatInfo[SatLabel][SatIdx["UDREI"]]) == 15:
                    SatCorrInfo["Flag"] = 0                             # Satellite Dont-Use (DU)

                # If the satellite can be used in PA service
                if SatCorrInfo["Flag"] == 1:
                
                    # CORRECTED ORBIT AND CLOCK
                    # -------------------------------------------------------
                    # Compute corrected Orbit and Clock, as well as the Sigma FLT

                    # Compute the DTR
                    SatOrb = np.array([float(SatInfo[SatLabel][SatIdx["SAT-X"]]),float(SatInfo[SatLabel][SatIdx["SAT-Y"]]),\
                        float(SatInfo[SatLabel][SatIdx["SAT-Z"]])])
                    SatVel = np.array([float(SatInfo[SatLabel][SatIdx["VEL-X"]]),float(SatInfo[SatLabel][SatIdx["VEL-Y"]]),\
                        float(SatInfo[SatLabel][SatIdx["VEL-Z"]])])
                    Dtr = - 2*(np.dot(SatOrb,SatVel)/Const.SPEED_OF_LIGHT**2)

                    # Corrected SBAS Orbit and Clock
                    SatCorrInfo["SatX"] = float(SatInfo[SatLabel][SatIdx["SAT-X"]]) + float(SatInfo[SatLabel][SatIdx["LTC-X"]])
                    SatCorrInfo["SatY"] = float(SatInfo[SatLabel][SatIdx["SAT-Y"]]) + float(SatInfo[SatLabel][SatIdx["LTC-Y"]])
                    SatCorrInfo["SatZ"] = float(SatInfo[SatLabel][SatIdx["SAT-Z"]]) + float(SatInfo[SatLabel][SatIdx["LTC-Z"]])
                    SatCorrInfo["SatClk"] = float(SatInfo[SatLabel][SatIdx["SAT-CLK"]]) + float(SatInfo[SatLabel][SatIdx["FC"]]) + \
                        float(SatInfo[SatLabel][SatIdx["LTC-B"]]) - float(SatInfo[SatLabel][SatIdx["TGD"]]) + Dtr

                    # Sigma FLT taking into account the degradation parameters
                    if int(SatInfo[SatLabel][SatIdx["RSS"]]) == 0:      
                        SatCorrInfo["SigmaFlt"] = (float(SatInfo[SatLabel][SatIdx["SIGMAUDRE"]])*float(SatInfo[SatLabel][SatIdx["DELTAUDRE"]]) + \
                            float(SatInfo[SatLabel][SatIdx["EPS-FC"]]) + float(SatInfo[SatLabel][SatIdx["EPS-RRC"]]) + \
                            float(SatInfo[SatLabel][SatIdx["EPS-LTC"]]) + float(SatInfo[SatLabel][SatIdx["EPS-ER"]]))**2
                    else:                                              
                        SatCorrInfo["SigmaFlt"] = (float(SatInfo[SatLabel][SatIdx["SIGMAUDRE"]])*float(SatInfo[SatLabel][SatIdx["DELTAUDRE"]]))**2 + \
                            float(SatInfo[SatLabel][SatIdx["EPS-FC"]])**2 + float(SatInfo[SatLabel][SatIdx["EPS-RRC"]])**2 + \
                            float(SatInfo[SatLabel][SatIdx["EPS-LTC"]])**2 + float(SatInfo[SatLabel][SatIdx["EPS-ER"]])**2
                
                    # CORRECTED TROPOSPHERIC DELAY
                    # -------------------------------------------------------
                    # Compute corrected Slant Tropospheric Delay, as well as the Sigma Tropo

                    # Compute the STD
                    SatCorrInfo["Std"] = LosInfo[SatLabel][LosIdx["STD"]]

                    # Compute the Sigma Tropo
                    Tpp = Tropo.computeTropoMappingFunction(SatCorrInfo["Elevation"])
                    SatCorrInfo["SigmaTropo"] = 0.12*Tpp


            # Prepare output for the satellite
            CorrInfo[SatLabel] = SatCorrInfo

        # End of if(SatPrepro["Status"] == 1):

    # End of for SatLabel, SatPrepro in PreproObsInfo.items():

    return CorrInfo
