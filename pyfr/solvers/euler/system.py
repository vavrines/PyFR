from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.euler.elements import EulerElements
from pyfr.solvers.euler.inters import (EulerIntInters, EulerPintInters,
                                       EulerMPIInters, EulerBaseBCInters)


class EulerSystem(BaseAdvectionSystem):
    name = 'euler'

    elementscls = EulerElements
    intinterscls = EulerIntInters
    pintinterscls = EulerPintInters
    mpiinterscls = EulerMPIInters
    bbcinterscls = EulerBaseBCInters
