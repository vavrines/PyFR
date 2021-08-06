# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.kpp.elements import KPPElements
from pyfr.solvers.kpp.inters import (KPPBaseBCInters,
                                           KPPIntInters,
                                           KPPMPIInters)


class KPPSystem(BaseAdvectionDiffusionSystem):
    name = 'kpp'

    elementscls = KPPElements
    intinterscls = KPPIntInters
    mpiinterscls = KPPMPIInters
    bbcinterscls = KPPBaseBCInters
