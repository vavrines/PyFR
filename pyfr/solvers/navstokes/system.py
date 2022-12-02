from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.navstokes.elements import NavierStokesElements
from pyfr.solvers.navstokes.inters import (NavierStokesBaseBCInters,
                                           NavierStokesIntInters,
                                           NavierStokesPintInters,
                                           NavierStokesMPIInters)


class NavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'navier-stokes'

    elementscls = NavierStokesElements
    intinterscls = NavierStokesIntInters
    pintinterscls = NavierStokesPintInters
    mpiinterscls = NavierStokesMPIInters
    bbcinterscls = NavierStokesBaseBCInters
