# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.bgk.elements import BGKElements
from pyfr.solvers.bgk.inters import (BGKIntInters, BGKMPIInters,
                                     BGKBaseBCInters)


class BGKSystem(BaseAdvectionSystem):
    name = 'bgk'

    elementscls = BGKElements
    intinterscls = BGKIntInters
    mpiinterscls = BGKMPIInters
    bbcinterscls = BGKBaseBCInters
