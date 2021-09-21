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
from math import sqrt, exp, cos, pi
import numpy as np

def IonoIGPs(Index, LosInfo):

    # Function defining the IGPs for the ionospheric interpolation

    NorEast = [float(LosInfo[LosIdx["IGP_NE_LON"]]), float(LosInfo[LosIdx["IGP_NE_LAT"]]), float(LosInfo[LosIdx["GIVD_NE"]]), float(LosInfo[LosIdx["GIVE_NE"]])]
    NorWest = [float(LosInfo[LosIdx["IGP_NW_LON"]]), float(LosInfo[LosIdx["IGP_NW_LAT"]]), float(LosInfo[LosIdx["GIVD_NW"]]), float(LosInfo[LosIdx["GIVE_NW"]])]
    SouWest = [float(LosInfo[LosIdx["IGP_SW_LON"]]), float(LosInfo[LosIdx["IGP_SW_LAT"]]), float(LosInfo[LosIdx["GIVD_SW"]]), float(LosInfo[LosIdx["GIVE_SW"]])]
    SouEast = [float(LosInfo[LosIdx["IGP_SE_LON"]]), float(LosInfo[LosIdx["IGP_SE_LAT"]]), float(LosInfo[LosIdx["GIVD_SE"]]), float(LosInfo[LosIdx["GIVE_SE"]])]
    IgpMask = [NorEast, NorWest, SouWest, SouEast]
    
    ActIgps = OrderedDict({})
    n = 0
    Counter = Index
    # Rectangular interpolation
    if Index == 0:
        while len(ActIgps.keys()) < 4:
            ActIgps[str(Counter + 1)] = IgpMask[Counter]
            Counter = Counter + 1
        return ActIgps 

    # Triangular interpolation
    if Index > 0: 
        while len(ActIgps.keys()) < 3:
            # Check if index has reached the end of the IgpMask list
            if Counter//(len(IgpMask)) == 1:
                Counter = 0
            ActIgps[str(n + 1)] = IgpMask[Counter]
            n = n + 1
            Counter = Counter + 1
        return ActIgps

def IonoInterpolation(LosInfo):

    # Function cpmputing the UIVD and the sigma UIRE for each line of sight (IPP) by interpolating the values 
    # of the surrounding active IGPs

    # Internal indexes definition
    InterIdx = {"LON": 0, "LAT": 1, "GIVD": 2, "GIVE": 3, "WEIGHT": 4}

    # Internal parameters definition
    Index = int(LosInfo[LosIdx["INTERP"]])
    IppLong = float(LosInfo[LosIdx["IPPLON"]])
    IppLat = float(LosInfo[LosIdx["IPPLAT"]])

    # Obtain active IGPs for the interpolation
    ActIgps = IonoIGPs(Index, LosInfo)

    # Compute the weights for the interpolation
    # Rectangular interpolation
    if Index == 0:
        if abs(IppLat) < 85.0:
            # IPPs between parallels N85 and S85         
            xpp = (IppLong - ActIgps["3"][InterIdx["LON"]])/(ActIgps["1"][InterIdx["LON"]] - ActIgps["3"][InterIdx["LON"]])
            ypp = (IppLat - ActIgps["3"][InterIdx["LAT"]])/(ActIgps["1"][InterIdx["LAT"]] - ActIgps["3"][InterIdx["LAT"]])
        else:
            # IPPs north of N85 or south of S85
            ypp = (abs(IppLat) - 85.0)/10
            xpp = ((IppLong - ActIgps["3"][InterIdx["LON"]])/90)*(1 - 2*ypp) + ypp
        # Append the weights to ActIgps dictionary
        ActIgps["1"].append(xpp*ypp)
        ActIgps["2"].append((1-xpp)*ypp)
        ActIgps["3"].append((1-xpp)*(1-ypp))
        ActIgps["4"].append(xpp*(1-ypp))

    # Triangular interpolation
    elif Index > 0: 
        if abs(IppLat) < 85.0:
            # IPPs between parallels N85 and S85
            if Index % 2 == 0:
                xpp = (IppLong - ActIgps["2"][InterIdx["LON"]])/(ActIgps["1"][InterIdx["LON"]] - ActIgps["2"][InterIdx["LON"]])
                ypp = (IppLat - ActIgps["2"][InterIdx["LAT"]])/(ActIgps["3"][InterIdx["LAT"]] - ActIgps["2"][InterIdx["LAT"]])
            else:
                xpp = (IppLong - ActIgps["2"][InterIdx["LON"]])/(ActIgps["3"][InterIdx["LON"]] - ActIgps["2"][InterIdx["LON"]])
                ypp = (IppLat - ActIgps["2"][InterIdx["LAT"]])/(ActIgps["1"][InterIdx["LAT"]] - ActIgps["2"][InterIdx["LAT"]])
        else:
            # IPPs north of N85 or south of S85
            LambdaWest = [ActIgps["1"][InterIdx["LON"]], ActIgps["3"][InterIdx["LON"]]]
            ypp = (abs(IppLat) - 85.0)/10
            xpp = ((IppLong - min(LambdaWest))/90)*(1 - 2*ypp) + ypp
        # Append the weights to ActIgps dictionary
        ActIgps["1"].append(ypp)
        ActIgps["2"].append(1-xpp-ypp)
        ActIgps["3"].append(xpp)

    # Compute UIVD at the IP
    UIVD = 0
    for i in range(1, len(ActIgps) + 1):
        UIVD = UIVD + ActIgps[str(i)][InterIdx["WEIGHT"]]*ActIgps[str(i)][InterIdx["GIVD"]]

    # Compute sigma UIRE
    UIVE2 = 0
    for i in range(1, len(ActIgps) + 1):
        UIVE2 = UIVE2 + ActIgps[str(i)][InterIdx["WEIGHT"]]*ActIgps[str(i)][InterIdx["GIVE"]]**2
    UIVE = sqrt(UIVE2)

    return UIVD, UIVE

def TropoInterpolation(RcvrInfo, Doy):

    # Function computing the STD for each line of sight by interpolating from MOPS values 

    # Internal indexes definition
    InterIdx = {"PPar": 0, "TPar": 1, "ePar": 2, "BetaPar": 3, "LambdaPar": 4}
    InterKey = {"Average": 0, "Variation": 1}

    # Meteorological data as a function of latitude
    MeteoParam = OrderedDict({})
    MeteoParam[15.0] = [[1013.25, 0.0], [299.65, 0.0], [26.31, 0.0], [6.3e-3, 0.0], [2.77, 0.0]]
    MeteoParam[30.0] = [[1017.25, -3.75], [294.15, 7.0], [21.79, 8.85], [6.05e-3, 0.25e-3], [3.15, 0.33]]
    MeteoParam[45.0] = [[1015.75, -2.25], [283.15, 11.0], [11.66, 7.24], [5.58e-3, 0.32e-3], [2.57, 0.46]]
    MeteoParam[60.0] = [[1011.75, -1.75], [272.15, 15.0], [6.78, 5.36], [5.39e-3, 0.81e-3], [1.81, 0.74]]
    MeteoParam[75.0] = [[1013.00, -0.50], [263.65, 14.5], [4.11, 3.39], [4.53e-3, 0.62e-3], [1.55, 0.3]]

    # Obtain local meteorological parameters around the receiver
    RcvrMeteoParam = []
    PrevLat = 15.0
    for Lat in MeteoParam.keys():
        if abs(float(RcvrInfo[RcvrIdx["LAT"]])) <= 15.0:
            RcvrMeteoParam = MeteoParam[15.0]
            break
        elif abs(float(RcvrInfo[RcvrIdx["LAT"]])) >= 75.0:
            RcvrMeteoParam = MeteoParam[75.0]
            break
        elif abs(float(RcvrInfo[RcvrIdx["LAT"]])) < Lat:
            for Idx in InterIdx.values():
                Parameters = []
                for Key in InterKey.values():
                    Parameters.append(MeteoParam[PrevLat][Idx][Key] + (MeteoParam[Lat][Idx][Key] - MeteoParam[PrevLat][Idx][Key]) * \
                        (abs((float(RcvrInfo[RcvrIdx["LAT"]])) - PrevLat)/(Lat-PrevLat)))
                RcvrMeteoParam.append(Parameters)
            break
        else:
            PrevLat = Lat
    
    # Compute meteorological parameters around the receiver
    if float(RcvrInfo[RcvrIdx["LAT"]]) > 0:
        # For northen latitudes
        Dmin = 28.0
    else:
        # For southern latitudes
        Dmin = 211.0

    RcvrP = RcvrMeteoParam[InterIdx["PPar"]][InterKey["Average"]] - RcvrMeteoParam[InterIdx["PPar"]][InterKey["Variation"]] * \
        cos(2*pi*(Doy - Dmin)/365.25)
    RcvrT = RcvrMeteoParam[InterIdx["TPar"]][InterKey["Average"]] - RcvrMeteoParam[InterIdx["TPar"]][InterKey["Variation"]] * \
        cos(2*pi*(Doy - Dmin)/365.25)
    Rcvre = RcvrMeteoParam[InterIdx["ePar"]][InterKey["Average"]] - RcvrMeteoParam[InterIdx["ePar"]][InterKey["Variation"]] * \
        cos(2*pi*(Doy - Dmin)/365.25)
    RcvrB = RcvrMeteoParam[InterIdx["BetaPar"]][InterKey["Average"]] - RcvrMeteoParam[InterIdx["BetaPar"]][InterKey["Variation"]] * \
        cos(2*pi*(Doy - Dmin)/365.25)
    RcvrL = RcvrMeteoParam[InterIdx["LambdaPar"]][InterKey["Average"]] - RcvrMeteoParam[InterIdx["LambdaPar"]][InterKey["Variation"]] * \
        cos(2*pi*(Doy - Dmin)/365.25)

    # Compute zero-altitude zenith delays
    Zhyd = 1e-6*((77.604*287.054)/9.784)*RcvrP
    Zwet = 1e-6*((382000*287.054)/(9.784*(RcvrL + 1) - RcvrB*287.054))*(Rcvre/RcvrT)

    # Compute dry and wet tropo delays
    Dhyd = Zhyd*(1-(RcvrB*float(RcvrInfo[RcvrIdx["ALT"]])/RcvrT))**(9.80665/(287.054*RcvrB))
    Dwet = Zwet*(1-(RcvrB*float(RcvrInfo[RcvrIdx["ALT"]])/RcvrT))**((((RcvrL + 1)*9.80665)/(287.054*RcvrB))-1)
    ZTD = Dhyd + Dwet

    return ZTD    

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
                "Flag": 0,              # 0: Not Used 1: Used for PA 2: Used for NPA
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
                    Dtr = - 2*(np.dot(SatOrb,SatVel)/Const.SPEED_OF_LIGHT)

                    # Corrected SBAS Orbit and Clock
                    SatCorrInfo["SatX"] = float(SatInfo[SatLabel][SatIdx["SAT-X"]]) + float(SatInfo[SatLabel][SatIdx["LTC-X"]])
                    SatCorrInfo["SatY"] = float(SatInfo[SatLabel][SatIdx["SAT-Y"]]) + float(SatInfo[SatLabel][SatIdx["LTC-Y"]])
                    SatCorrInfo["SatZ"] = float(SatInfo[SatLabel][SatIdx["SAT-Z"]]) + float(SatInfo[SatLabel][SatIdx["LTC-Z"]])
                    SatCorrInfo["SatClk"] = float(SatInfo[SatLabel][SatIdx["SAT-CLK"]]) + float(SatInfo[SatLabel][SatIdx["FC"]]) + \
                        float(SatInfo[SatLabel][SatIdx["LTC-B"]]) - float(SatInfo[SatLabel][SatIdx["TGD"]]) + Dtr

                    # Sigma FLT taking into account the degradation parameters
                    if int(SatInfo[SatLabel][SatIdx["RSS"]]) == 0:      
                        SatCorrInfo["SigmaFlt"] = float(SatInfo[SatLabel][SatIdx["SIGMAUDRE"]])*float(SatInfo[SatLabel][SatIdx["DELTAUDRE"]]) + \
                            float(SatInfo[SatLabel][SatIdx["EPS-FC"]]) + float(SatInfo[SatLabel][SatIdx["EPS-RRC"]]) + \
                            float(SatInfo[SatLabel][SatIdx["EPS-LTC"]]) + float(SatInfo[SatLabel][SatIdx["EPS-ER"]])
                    else:                                              
                        SatCorrInfo["SigmaFlt"] = sqrt((float(SatInfo[SatLabel][SatIdx["SIGMAUDRE"]])*float(SatInfo[SatLabel][SatIdx["DELTAUDRE"]]))**2 + \
                            float(SatInfo[SatLabel][SatIdx["EPS-FC"]])**2 + float(SatInfo[SatLabel][SatIdx["EPS-RRC"]])**2 + \
                            float(SatInfo[SatLabel][SatIdx["EPS-LTC"]])**2 + float(SatInfo[SatLabel][SatIdx["EPS-ER"]])**2)
                
                    # IONOSPHERIC DELAY
                    # -------------------------------------------------------
                    # Compute corrected Slant Ionospheric Delay UISD, as well as the Sigma UIRE

                    # Compute the UIVD and the Mpp
                    UIVD, UIVE = IonoInterpolation(LosInfo[SatLabel])
                    Fpp = Iono.computeIonoMappingFunction(SatCorrInfo["Elevation"])
                    
                    # Compute the UISD
                    SatCorrInfo["Uisd"] = Fpp*UIVD

                    # Compute the Sigma UIRE
                    SatCorrInfo["SigmaUire"] = Fpp*UIVE
                    
                    # TROPOSPHERIC DELAY
                    # -------------------------------------------------------
                    # Compute corrected Slant Tropospheric Delay STD, as well as the Sigma TROPO according to MOPS

                    # Compute the STD
                    Tpp = Tropo.computeTropoMappingFunction(SatCorrInfo["Elevation"])
                    SatCorrInfo["Std"] = TropoInterpolation(Rcvr, SatCorrInfo["Doy"])*Tpp
                    # SatCorrInfo["Std"] = float(LosInfo[SatLabel][LosIdx["STD"]])

                    # Compute the Sigma Tropo
                    SatCorrInfo["SigmaTropo"] = 0.12*Tpp

                    # USER AIRBORNE SIGMA
                    # -------------------------------------------------------
                    # Compute user Sigma AIR according to MOPS

                    # Compute Sigma Multipath
                    if SatCorrInfo["Elevation"] < 2:
                        # Display error
                        sys.stderr.write("ERROR: Elevation angle is lower than 2 degrees")
                        sys.exit(-1)
                    else:
                        SatCorrInfo["SigmaMultipath"] = 0.13 + 0.53*exp(-SatCorrInfo["Elevation"]/10)

                    # Compute Sigma Noise + Divergence 
                    if SatCorrInfo["Elevation"] < float(Conf["ELEV_NOISE_TH"]):
                        SatCorrInfo["SigmaNoiseDiv"] = 0.36
                    else:
                        SatCorrInfo["SigmaNoiseDiv"] = 0.15

                    # Compute Sigma Airborne
                    if int(Conf["EQUIPMENT_CLASS"]) == 2 or int(Conf["EQUIPMENT_CLASS"]) == 3 or int(Conf["EQUIPMENT_CLASS"]) == 4:
                        SatCorrInfo["SigmaAirborne"] = sqrt(SatCorrInfo["SigmaMultipath"]**2 + SatCorrInfo["SigmaNoiseDiv"]**2)

                    # USER SIGMA
                    # -------------------------------------------------------
                    # Compute user sigma UERE

                    SatCorrInfo["SigmaUere"] = sqrt(SatCorrInfo["SigmaFlt"]**2 + SatCorrInfo["SigmaUire"]**2 + SatCorrInfo["SigmaTropo"]**2 + \
                        SatCorrInfo["SigmaAirborne"]**2)

                    # CORRECTED PSEUDO-RANGE AND GEOMETRICAL RANGE
                    # -------------------------------------------------------
                    # Compute the corrected pseudo-range and the satellite geometrical range 
                    # Compute the first residuals

                    # Calculate corrected pseudo-range
                    SatCorrInfo["CorrPsr"] = SatPrepro["SmoothC1"] + SatCorrInfo["SatClk"] -  SatCorrInfo["Uisd"] - SatCorrInfo["Std"]

                    # Compute satellite geometrical range
                    RcvrPos = np.array(Rcvr[RcvrIdx["XYZ"]])
                    SatPos = np.array([SatCorrInfo["SatX"], SatCorrInfo["SatY"], SatCorrInfo["SatZ"]])
                    Sat_Rcvr = np.subtract(SatPos, RcvrPos)
                    SatCorrInfo["GeomRange"] = np.linalg.norm(Sat_Rcvr)

                    # Obtain first residuals
                    SatCorrInfo["PsrResidual"] = SatCorrInfo["CorrPsr"] - SatCorrInfo["GeomRange"]  

                    # UPDATE SOME PARAMETERS
                    # -------------------------------------------------------      

                    # Satellite LOS unitary vector and instantaneous orbit corrections vector
                    ULos = Sat_Rcvr/SatCorrInfo["GeomRange"]
                    LTCxyz = np.array([float(SatInfo[SatLabel][SatIdx["LTC-X"]]), float(SatInfo[SatLabel][SatIdx["LTC-Y"]]), float(SatInfo[SatLabel][SatIdx["LTC-Z"]])])

                    # Weighted sum of the residuals
                    ResN = ResN + 1/(SatCorrInfo["SigmaUere"])**2
                    ResSum = ResSum + SatCorrInfo["PsrResidual"]/(SatCorrInfo["SigmaUere"])**2

                    # ENT-GPS offset parameters
                    EntGpsN = EntGpsN + 1
                    EntGpsSum = EntGpsSum + (np.dot(LTCxyz, ULos) - (float(SatInfo[SatLabel][SatIdx["FC"]]) + float(SatInfo[SatLabel][SatIdx["LTC-B"]])))

                # End of if SatCorrInfo["Flag"] == 1:

            # Prepare output for the satellite
            CorrInfo[SatLabel] = SatCorrInfo

        # End of if SatPrepro["Status"] == 1:

    # End of for SatLabel, SatPrepro in PreproObsInfo.items():

    # FINAL RESIDUALS AND ENT-GPS OFFSET
    # -------------------------------------------------------
    # Compute the final residuals by substracting the guess of the receiver clock
    # Compute the ENT to GPS time offset

    # Check that there are at least one monitored satellite at current epoch
    if EntGpsN > 0:
        for SatLabel, SatCorrInfo in CorrInfo.items():
            # Compute ENT-GPS offset
            SatCorrInfo["EntGps"] = EntGpsSum/EntGpsN
            # Compute the final residuals
            SatCorrInfo["RcvrClk"] = ResSum/ResN
            # Compute the final residuals
            SatCorrInfo["PsrResidual"] = SatCorrInfo["PsrResidual"] - SatCorrInfo["RcvrClk"]

    return CorrInfo
