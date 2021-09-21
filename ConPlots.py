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
ConfPlots["PLOT_SAT_TRACKS"] = 1
ConfPlots["PLOT_LTC"] = 1
ConfPlots["PLOT_ENT_GPS"] = 1
ConfPlots["PLOT_SIGMA_FLT"] = 1
ConfPlots["PLOT_UIVD"] = 1
ConfPlots["PLOT_UISD"] = 1
ConfPlots["PLOT_SIGMA_UIRE"] = 1
ConfPlots["PLOT_STD"] = 1
ConfPlots["PLOT_SIGMA_TROPO"] = 1
ConfPlots["PLOT_SIGMA_MULTI"] = 1
ConfPlots["PLOT_SIGMA_NOISE"] = 1
ConfPlots["PLOT_SIGMA_AIR"] = 1
ConfPlots["PLOT_SIGMA_UERE"] = 1
ConfPlots["PLOT_SIGMA_UERE_STAT"] = 1
ConfPlots["PLOT_RCVR_CLK"] = 1
ConfPlots["PLOT_RES"] = 1
ConfPlots["PLOT_SIGMA_UERE_STATS"] = 1

