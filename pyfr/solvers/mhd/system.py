from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.mhd.elements import MHDElements
from pyfr.solvers.mhd.inters import (MHDIntInters, MHDMPIInters,
                                     MHDBaseBCInters)


class MHDSystem(BaseAdvectionSystem):
    name = 'mhd'

    elementscls = MHDElements
    intinterscls = MHDIntInters
    mpiinterscls = MHDMPIInters
    bbcinterscls = MHDBaseBCInters
