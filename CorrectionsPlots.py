#!/usr/bin/env python

########################################################################
# PETRUS/SRC/PreprocessingPlots.py:
# This is the PreprocessingPlots Module of PETRUS tool
#
#  Project:        PETRUS
#  File:           PreprocessingPlots.py
#  Date(YY/MM/DD): 05/02/21
#
#  Author: GNSS Academy
#  Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################

import sys, os
from pandas import unique
from pandas import read_csv
from InputOutput import CorrIdx, SatIdx
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants, Iono
from COMMON.Plots import generatePlot
from COMMON.Coordinates import xyz2llh
import numpy as np
from ConPlots import ConfPlots
import matplotlib.pyplot as plt

def initPlot(CorrFile, PlotConf, Title, Label):
    
    # Compute information from PreproObsFile
    CorrFileName = os.path.basename(CorrFile)
    CorrFileNameSplit = CorrFileName.split('_')
    Rcvr = CorrFileNameSplit[1]
    DatepDat = CorrFileNameSplit[2]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]

    # Dump information into PlotConf
    PlotConf["xLabel"] = "Hour of Day %s" % Doy

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/CORR/Figures/%s/' % Label + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)

# Plot Monitored Satellite Tracks
def plotSatTracks(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "Monitored Satellites Tracks", "SATS_TRACKS_vs_TIME")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (16.8,15.2)

    PlotConf["LonMin"] = -115
    PlotConf["LonMax"] = 125
    PlotConf["LatMin"] = -35
    PlotConf["LatMax"] = 65
    PlotConf["LonStep"] = 15
    PlotConf["LatStep"] = 10

    PlotConf["yTicks"] = range(PlotConf["LatMin"],PlotConf["LatMax"]+1,10)
    PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

    PlotConf["xTicks"] = range(PlotConf["LonMin"],PlotConf["LonMax"]+1,15)
    PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Processing data to be plotted
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    x = CorrData[CorrIdx["SAT-X"]][FilterCond].to_numpy()
    y = CorrData[CorrIdx["SAT-Y"]][FilterCond].to_numpy()
    z = CorrData[CorrIdx["SAT-Z"]][FilterCond].to_numpy()
    Longitude = np.zeros(len(x))
    Latitude = np.zeros(len(x))

    for index in range(len(x)):
        Longitude[index], Latitude[index], h = xyz2llh(x[index], y[index], z[index])

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [Deg]"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 90.0
    PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0
    PlotConf["xData"][Label] = Longitude
    PlotConf["yData"][Label] = Latitude
    PlotConf["zData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot LTC Corrections
def plotLtcCorr(CorrFile, CorrData, SatData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "LTC Corrections", "LTC_CORRECTIONS_vs_TIME")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Background"] = ["LTC-X","LTC-Y","LTC-Z","LTC-B","FC"]

    PlotConf["yLabel"] = "Fast and Long Term Corrections [m]"

    PlotConf["xTicks"] = range(0,25)
    PlotConf["xLim"] = [0,24]

    PlotConf["Grid"] = True
    PlotConf["Legend"] = ["ENTtoGPS","LTC-X","LTC-Y","LTC-Z","LTC-B","FC"]

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Plotting
    PlotConf["xData"], PlotConf["xData2"] = {}, {}
    PlotConf["yData"], PlotConf["yData2"] = {}, {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    for Label in PlotConf["Legend"]:
        if Label == "ENTtoGPS":
            PlotConf["Marker"] = '_'
            PlotConf["LineWidth"] = 1.0
            PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
            PlotConf["yData"][Label] = CorrData[CorrIdx[Label]][FilterCond]
        else:
            PlotConf["xData"][Label] = SatData[SatIdx["SOD"]] / GnssConstants.S_IN_H
            PlotConf["yData"][Label] = SatData[SatIdx[Label]]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot ENT-GPS Offset
def plotEntGps(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "ENT - GPS Offset", "ENT_GPS_OFFSET_vs_TIME")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)

    PlotConf["yLabel"] = "ENT-GPS Offset [m]"

    PlotConf["xTicks"] = range(0,25)
    PlotConf["xLim"] = [0,24]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '+'
    PlotConf["LineWidth"] = 0.75

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["SOD"]][FilterCond] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = CorrData[CorrIdx["ENTtoGPS"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Sigma FLT
def plotSigmaFlt(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "Sigma FLT", "SIGMA_FLT_vs_ELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)

    PlotConf["yLabel"] = "Sigma FLT [m]"

    PlotConf["xLabel"] = "Elevation [Deg]"
    PlotConf["xTicks"] = range(0,100,10)
    PlotConf["xLim"] = [0,90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PrnList = sorted(unique(CorrData[CorrIdx["PRN"]]))
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = min(PrnList)
    PlotConf["ColorBarMax"] = max(PrnList)
    PlotConf["ColorBarTicks"] = range(min(PrnList), max(PrnList) + 1)

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SFLT"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot UIVD
def plotUivd(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "UIVD", "UIVD_vs_TIME")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)

    PlotConf["yLabel"] = "UIVD [m]"

    PlotConf["xTicks"] = range(0,25)
    PlotConf["xLim"] = [0,24]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Processing data to be plotted
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Sod = CorrData[CorrIdx["SOD"]][FilterCond].to_numpy()
    Elev = CorrData[CorrIdx["ELEV"]][FilterCond].to_numpy()
    Uisd = CorrData[CorrIdx["UISD"]][FilterCond].to_numpy()
    Uivd = np.zeros(len(Sod))
    
    for index in range(len(Uivd)):
        Uivd[index] = Uisd[index]/Iono.computeIonoMappingFunction(Elev[index])

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [Deg]"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 90.0
    PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0
    PlotConf["xData"][Label] = Sod / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = Uivd
    PlotConf["zData"][Label] = Elev

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot UISD
def plotUisd(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "UISD", "UISD_at_IPP")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (16.8,15.2)

    PlotConf["LonMin"] = -20
    PlotConf["LonMax"] = 25
    PlotConf["LatMin"] = 25
    PlotConf["LatMax"] = 60
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["yTicks"] = range(PlotConf["LatMin"],PlotConf["LatMax"]+1,10)
    PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

    PlotConf["xTicks"] = range(PlotConf["LonMin"],PlotConf["LonMax"]+1,15)
    PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "UISD [m]"
    PlotConf["ColorBarMax"] = int(max(CorrData[CorrIdx["UISD"]])) + 1
    PlotConf["ColorBarMin"] = int(min(CorrData[CorrIdx["UISD"]]))
    PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["IPPLON"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["IPPLAT"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["UISD"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Sigma UIRE
def plotSigmaUire(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "Sigma UIRE", "SIGMA_UIRE_vs_ELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)

    PlotConf["yLabel"] = "Sigma UIRE [m]"

    PlotConf["xLabel"] = "Elevation [Deg]"
    PlotConf["xTicks"] = range(0,100,10)
    PlotConf["xLim"] = [0,90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PrnList = sorted(unique(CorrData[CorrIdx["PRN"]]))
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = min(PrnList)
    PlotConf["ColorBarMax"] = max(PrnList)
    PlotConf["ColorBarTicks"] = range(min(PrnList), max(PrnList) + 1)

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SUIRE"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot STD
def plotStd(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "STD", "STD_at_IPP")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (16.8,15.2)

    PlotConf["LonMin"] = -20
    PlotConf["LonMax"] = 25
    PlotConf["LatMin"] = 25
    PlotConf["LatMax"] = 60
    PlotConf["LonStep"] = 5
    PlotConf["LatStep"] = 5

    PlotConf["yTicks"] = range(PlotConf["LatMin"],PlotConf["LatMax"]+1,10)
    PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

    PlotConf["xTicks"] = range(PlotConf["LonMin"],PlotConf["LonMax"]+1,15)
    PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

    PlotConf["Grid"] = True
    PlotConf["Map"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "STD [m]"
    PlotConf["ColorBarMax"] = int(max(CorrData[CorrIdx["STD"]])) + 1
    PlotConf["ColorBarMin"] = int(min(CorrData[CorrIdx["STD"]]))
    PlotConf["ColorBarTicks"] = None

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["IPPLON"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["IPPLAT"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["STD"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Sigma TROPO
def plotSigmaTropo(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "Sigma TROPO", "SIGMA_TROPO_vs_ELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)

    PlotConf["yLabel"] = "Sigma TROPO [m]"

    PlotConf["xLabel"] = "Elevation [Deg]"
    PlotConf["xTicks"] = range(0,100,10)
    PlotConf["xLim"] = [0,90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PrnList = sorted(unique(CorrData[CorrIdx["PRN"]]))
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = min(PrnList)
    PlotConf["ColorBarMax"] = max(PrnList)
    PlotConf["ColorBarTicks"] = range(min(PrnList), max(PrnList) + 1)

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["STROPO"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Sigma Multipath
def plotSigmaMulti(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "Sigma Multipath", "SIGMA_MULTI_vs_ELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)

    PlotConf["yLabel"] = "Sigma Multipath [m]"

    PlotConf["xLabel"] = "Elevation [Deg]"
    PlotConf["xTicks"] = range(0,100,10)
    PlotConf["xLim"] = [0,90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PrnList = sorted(unique(CorrData[CorrIdx["PRN"]]))
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = min(PrnList)
    PlotConf["ColorBarMax"] = max(PrnList)
    PlotConf["ColorBarTicks"] = range(min(PrnList), max(PrnList) + 1)

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SMP"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Sigma Noise + Divergence
def plotSigmaNoise(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "Sigma Noise + Divergence", "SIGMA_NOISE_vs_ELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)

    PlotConf["yLabel"] = "Sigma Noise + Divergence [m]"

    PlotConf["xLabel"] = "Elevation [Deg]"
    PlotConf["xTicks"] = range(0,100,10)
    PlotConf["xLim"] = [0,90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PrnList = sorted(unique(CorrData[CorrIdx["PRN"]]))
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = min(PrnList)
    PlotConf["ColorBarMax"] = max(PrnList)
    PlotConf["ColorBarTicks"] = range(min(PrnList), max(PrnList) + 1)

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SNOISEDIV"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Sigma Airborne
def plotSigmaAir(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "Sigma Airborne", "SIGMA_AIR_vs_ELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)

    PlotConf["yLabel"] = "Sigma Airborne [m]"

    PlotConf["xLabel"] = "Elevation [Deg]"
    PlotConf["xTicks"] = range(0,100,10)
    PlotConf["xLim"] = [0,90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PrnList = sorted(unique(CorrData[CorrIdx["PRN"]]))
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = min(PrnList)
    PlotConf["ColorBarMax"] = max(PrnList)
    PlotConf["ColorBarTicks"] = range(min(PrnList), max(PrnList) + 1)

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SAIR"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Sigma Airborne
def plotSigmaUere(CorrFile, CorrData):

    # Graph settings definition
    PlotConf = {}
    initPlot(CorrFile, PlotConf, "Sigma UERE", "SIGMA_UERE_vs_ELEV")

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)

    PlotConf["yLabel"] = "Sigma UERE [m]"

    PlotConf["xLabel"] = "Elevation [Deg]"
    PlotConf["xTicks"] = range(0,100,10)
    PlotConf["xLim"] = [0,90]

    PlotConf["Grid"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.75

    # Colorbar definition
    PrnList = sorted(unique(CorrData[CorrIdx["PRN"]]))
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = min(PrnList)
    PlotConf["ColorBarMax"] = max(PrnList)
    PlotConf["ColorBarTicks"] = range(min(PrnList), max(PrnList) + 1)

    # Plotting
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    FilterCond = CorrData[CorrIdx["FLAG"]] == 1
    Label = 0
    PlotConf["xData"][Label] = CorrData[CorrIdx["ELEV"]][FilterCond]
    PlotConf["yData"][Label] = CorrData[CorrIdx["SUERE"]][FilterCond]
    PlotConf["zData"][Label] = CorrData[CorrIdx["PRN"]][FilterCond]

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

def generateCorrPlots(CorrFile, SatFile, RcvrInfo):
    
    # Purpose: generate output plots regarding Preprocessing results

    # Parameters
    # ==========
    # CorrFile: str
    #           Path to CORR output file
    # SatFile: str
    #          Path to SAT input file
    # RcvrInfo: List
    #           List containing information on the receiver

    # Returns
    # =======
    # Nothing
    
    # Monitored Satellite Tracks
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_SAT_TRACKS"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SAT-X"],CorrIdx["SAT-Y"],CorrIdx["SAT-Z"],CorrIdx["ELEV"],CorrIdx["FLAG"]])

        print( 'Plot Monitored Satellites Tracks vs Time ...')
      
        # Configure plot and call plot generation function
        plotSatTracks(CorrFile, CorrData)

    # LTC Corrections
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_LTC"] == 1):
        # Read the cols we need from SatFile file
        SatData = read_csv(SatFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[SatIdx["SOD"],SatIdx["LTC-X"],SatIdx["LTC-Y"],SatIdx["LTC-Z"],SatIdx["LTC-B"],SatIdx["FC"]])

        # Read the cols we need from CorrFile file
        CorrData = read_csv(SatFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"],CorrIdx["ENTtoGPS"],CorrIdx["FLAG"]])

        print( 'Plot LTC corrections vs Time ...')
      
        # Configure plot and call plot generation function
        plotLtcCorr(CorrFile, CorrData, SatData)

    # ENT-GPS Offset
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_ENT_GPS"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"],CorrIdx["ENTtoGPS"],CorrIdx["FLAG"]])

        print( 'Plot ENT-GPS Offset vs Time ...')
      
        # Configure plot and call plot generation function
        plotEntGps(CorrFile, CorrData)

    # Sigma FLT
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_SIGMA_FLT"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["PRN"],CorrIdx["ELEV"],CorrIdx["SFLT"],CorrIdx["FLAG"]])

        print( 'Plot Sigma FLT vs Elevation ...')
      
        # Configure plot and call plot generation function
        plotSigmaFlt(CorrFile, CorrData)

    # UIVD
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_UIVD"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["SOD"],CorrIdx["ELEV"],CorrIdx["UISD"],CorrIdx["FLAG"]])

        print( 'Plot UIVD vs Time ...')
      
        # Configure plot and call plot generation function
        plotUivd(CorrFile, CorrData)

    # UISD
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_UISD"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["IPPLON"],CorrIdx["IPPLAT"],CorrIdx["UISD"],CorrIdx["FLAG"]])

        print( 'Plot UISD vs Time ...')
      
        # Configure plot and call plot generation function
        plotUisd(CorrFile, CorrData)

    # Sigma UIRE
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_SIGMA_UIRE"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["PRN"],CorrIdx["ELEV"],CorrIdx["SUIRE"],CorrIdx["FLAG"]])

        print( 'Plot Sigma UIRE vs Elevation ...')
      
        # Configure plot and call plot generation function
        plotSigmaUire(CorrFile, CorrData)

    # STD
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_STD"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["IPPLON"],CorrIdx["IPPLAT"],CorrIdx["STD"],CorrIdx["FLAG"]])

        print( 'Plot STD vs Time ...')
      
        # Configure plot and call plot generation function
        plotStd(CorrFile, CorrData)
    
    # Sigma TROPO
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_SIGMA_TROPO"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["PRN"],CorrIdx["ELEV"],CorrIdx["STROPO"],CorrIdx["FLAG"]])

        print( 'Plot Sigma TROPO vs Elevation ...')
      
        # Configure plot and call plot generation function
        plotSigmaTropo(CorrFile, CorrData)
    
    # Sigma Multipath
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_SIGMA_MULTI"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["PRN"],CorrIdx["ELEV"],CorrIdx["SMP"],CorrIdx["FLAG"]])

        print( 'Plot Sigma Multipath vs Elevation ...')
      
        # Configure plot and call plot generation function
        plotSigmaMulti(CorrFile, CorrData)
    
    # Sigma Noise + Divergence
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_SIGMA_NOISE"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["PRN"],CorrIdx["ELEV"],CorrIdx["SNOISEDIV"],CorrIdx["FLAG"]])

        print( 'Plot Sigma Noise + Divergence vs Elevation ...')
      
        # Configure plot and call plot generation function
        plotSigmaNoise(CorrFile, CorrData)
    
    # Sigma AIR
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_SIGMA_AIR"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["PRN"],CorrIdx["ELEV"],CorrIdx["SAIR"],CorrIdx["FLAG"]])

        print( 'Plot Sigma Airborne vs Elevation ...')
      
        # Configure plot and call plot generation function
        plotSigmaAir(CorrFile, CorrData)
    
    # Sigma UERE
    # ----------------------------------------------------------
    if(ConfPlots["PLOT_SIGMA_UERE"] == 1):
        # Read the cols we need from CorrFile file
        CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
        usecols=[CorrIdx["PRN"],CorrIdx["ELEV"],CorrIdx["SUERE"],CorrIdx["FLAG"]])

        print( 'Plot Sigma UERE vs Elevation ...')
      
        # Configure plot and call plot generation function
        plotSigmaUere(CorrFile, CorrData)