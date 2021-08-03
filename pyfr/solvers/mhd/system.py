# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.mhd.elements import MHDElements
from pyfr.solvers.mhd.inters import (MHDBaseBCInters,
                                           MHDIntInters,
                                           MHDMPIInters)


class MHDSystem(BaseAdvectionDiffusionSystem):
    name = 'mhd'

    elementscls = MHDElements
    intinterscls = MHDIntInters
    mpiinterscls = MHDMPIInters
    bbcinterscls = MHDBaseBCInters
