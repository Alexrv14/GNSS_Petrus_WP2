#!/usr/bin/env python

########################################################################
# PETRUS/SRC/Preprocessing.py:
# This is the Inputs (conf and input files) Module of PETRUS tool
#
#  Project:        PETRUS
#  File:           Preprocessing.py
#  Date(YY/MM/DD): 01/02/21
#
#  Author: GNSS Academy
#  Copyright 2021 GNSS Academy
#
# -----------------------------------------------------------------
# Date       | Author             | Action
# -----------------------------------------------------------------
#
########################################################################

# Import External and Internal functions and Libraries
#----------------------------------------------------------------------
from collections import OrderedDict

# Plots configuration flags
ConfPlots = OrderedDict({})
ConfPlots["PLOT_SAT_TRACKS"] = 0
ConfPlots["PLOT_LTC"] = 0
ConfPlots["PLOT_ENT_GPS"] = 0
ConfPlots["PLOT_SIGMA_FLT"] = 0
ConfPlots["PLOT_UIVD"] = 0
ConfPlots["PLOT_UISD"] = 0
ConfPlots["PLOT_SIGMA_UIRE"] = 0
ConfPlots["PLOT_STD"] = 0
ConfPlots["PLOT_SIGMA_TROPO"] = 0
ConfPlots["PLOT_SIGMA_MULTI"] = 0
ConfPlots["PLOT_SIGMA_NOISE"] = 0
ConfPlots["PLOT_SIGMA_AIR"] = 0
ConfPlots["PLOT_SIGMA_UERE"] = 0
ConfPlots["PLOT_SIGMA_UERE_STAT"] = 0
ConfPlots["PLOT_RCVR_CLK"] = 0
ConfPlots["PLOT_RES"] = 0
ConfPlots["PLOT_SIGMA_UERE_STATS"] = 0

