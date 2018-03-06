import math
import numpy as np
from Input_variables import *
'''
It finds the position and the area of the booms
'''
# in cm
NSpS = (n_st - 1) / 2
CirC = h_a * math.pi
Lsl = math.sqrt((c_a - h_a / 2) ** 2 + (h_a / 2) ** 2)
# h_alf of the surface distance of the Aileron
Sa = CirC / 4 + Lsl
# spacing of stringers
StrS = Sa / (NSpS + 1)
# Stringer loc_ation along Surface
StArea = h_st * t_st + (w_st - t_st) * t_st
SpA = h_a * t_sp
Ac = (h_a / 2) ** 2 * math.pi
At = h_a * (c_a - h_a / 2) / 2
Ax = At + Ac

# Surface position of stringers
# need if statement to determine if on the curve or the slope
def StringerSurface(NSpS, StrS):
    StrinLoc = np.zeros(NSpS + 2)
    for i in xrange(0, NSpS, 1):  # Stringer Loc_ation along the Surface
        StrinLoc[i + 1] = StrS * (i + 1)
    return StrinLoc
    # Stringer loc_ation Z-Y plane


StrinLoc = StringerSurface(NSpS, StrS)


def StringerZ(h_a, c_a, NSpS, CirC, Lsl, StrinLoc):
    StZ = np.zeros(NSpS + 2)
    for z in xrange(0, NSpS, 1):  # stringer loc_ation on the Z plane
        if StrinLoc[z] < CirC / 4:  # In curve
            StZ[z + 1] = h_a / 2 * (1 - math.cos(StrinLoc[z + 1] / (CirC / 4) * math.pi / 2))
        elif StrinLoc[z] >= CirC / 4:  # in slope
            StZ[z + 1] = (StrinLoc[z + 1] - CirC / 4) * (c_a - h_a / 2) / Lsl + h_a / 2
    return StZ


StZ = StringerZ(h_a, c_a, NSpS, CirC, Lsl, StrinLoc)


def StringerY(h_a, NSpS, CirC, Lsl, StrinLoc):
    Sty = np.zeros(NSpS + 2)
    for y in xrange(0, NSpS, 1):  # stringer loc_ation on the Y plane
        if StrinLoc[y] < CirC / 4:  # in circle
            Sty[y + 1] = h_a / 2 * (math.sin(StrinLoc[y + 1] / (CirC / 4) * math.pi / 2))
        elif StrinLoc[y] >= CirC / 4:  # in slope
            Sty[y + 1] = (Lsl - (StrinLoc[y + 1] - CirC / 4)) * (h_a / 2) / Lsl
    return Sty


Sty = StringerY(h_a, NSpS, CirC, Lsl, StrinLoc)


# in cm
# area of spar
# Boom Area###########3
# Due to stringer Area, including an extra boom due to the spar running through the cross-section

# define the boom loc_ation for the spar at h_a/2 for Z and Y


# determining the Boom area due to the skin thickness.
# due to assumption th_at the normal forces on the top part and the bottom section are
# equal but opposite, then the area for the top boom and its mirrored counter part are the same.
# the areas of the booms th_at are loc_ated at different z and y loc_ations are determined by their distance
# to the centroid, based on their ratio, so sigma1/sigma2 would just be y1/y2
# 9 boom loc_ations per side, 8 h_ave a skin contribution due to the initial stringer lying on the neutral axis.

def BoomArea(h_a, NSpS, CirC, StArea, t_sk, t_sp, Sty, StZ,
             StrinLoc):  # the boom area for each boom as specified in the report.
    BrA = np.zeros(NSpS + 2)
    StrinLoc[NSpS + 1] = CirC / 4.0
    StZ[NSpS + 1] = h_a / 2.0
    Sty[NSpS + 1] = h_a / 2.0

    BrA[0] = StArea

    BrA[1] = StArea + t_sk * (StrinLoc[1] - StrinLoc[0]) / 6 * (2)

    BrA[8] = t_sp * (h_a) / 6 * (1) + t_sk * (StrinLoc[8] - StrinLoc[1]) / 6 * (2 + Sty[1] / Sty[8]) + t_sk * (
            StrinLoc[2] - StrinLoc[8]) / 6 * (2 + Sty[2] / Sty[8])

    BrA[2] = StArea + t_sk * (StrinLoc[2] - StrinLoc[8]) / 6 * (2 + Sty[8] / Sty[2]) + t_sk * (
            StrinLoc[3] - StrinLoc[2]) / 6 * (2 + Sty[3] / Sty[2])

    BrA[3] = StArea + t_sk * (StrinLoc[3] - StrinLoc[2]) / 6 * (2 + Sty[2] / Sty[3]) + t_sk * (
            StrinLoc[4] - StrinLoc[3]) / 6 * (2 + Sty[3] / Sty[4])

    BrA[4] = StArea + t_sk * (StrinLoc[4] - StrinLoc[3]) / 6 * (2 + Sty[3] / Sty[4]) + t_sk * (
            StrinLoc[5] - StrinLoc[4]) / 6 * (2 + Sty[4] / Sty[5])

    BrA[5] = StArea + t_sk * (StrinLoc[5] - StrinLoc[4]) / 6 * (2 + Sty[4] / Sty[5]) + t_sk * (
            StrinLoc[6] - StrinLoc[5]) / 6 * (2 + Sty[5] / Sty[6])

    BrA[6] = StArea + t_sk * (StrinLoc[6] - StrinLoc[5]) / 6 * (2 + Sty[5] / Sty[6]) + t_sk * (
            StrinLoc[7] - StrinLoc[6]) / 6 * (2 + Sty[6] / Sty[7])

    BrA[7] = StArea + t_sk * (StrinLoc[7] - StrinLoc[6]) / 6 * (2 + Sty[6] / Sty[7]) + t_sk * 2 * (
            Sa - StrinLoc[7]) / 6 * (2)

    return BrA


BrA = BoomArea(h_a, NSpS, CirC, StArea, t_sk, t_sp, Sty, StZ, StrinLoc)
