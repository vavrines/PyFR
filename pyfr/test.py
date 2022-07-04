from pyfr.shapes import *
from pyfr.polys import *
from pyfr.inifile import Inifile

from shapes import TetLShape



path = '/mnt/c/Users/Tarik-Personal/Documents/GitHub/PyFR/pyfr/test.ini'
cfg = Inifile.load(path)


def test_m0():

    pyr = HexShape(False, cfg)

    m0 = pyr.m0
    tol = 1e-6

    # Test 1: Identity mats for collocated points is recovered
    for i in range(pyr.nupts):
        row = m0[i,:]
        assert np.abs(np.max(row) - 1) < tol
        assert np.abs(np.min(row)) < tol
        assert np.abs(np.sum(np.abs(row)) - 1) < tol
    print('Passed M0 test 1: Identity')

    # Test 2: Recover constant
    u = np.ones(pyr.nupts)
    assert np.amax(m0 @ u - np.ones(pyr.nfpts)) < tol
    print('Passed M0 test 2: Constant states')

    # Test 2: Recover linear
    x = pyr.upts[:,0]
    y = pyr.upts[:,1]
    z = pyr.upts[:,2]
    assert np.amax(m0 @ x - pyr.fpts[:,0]) < tol
    assert np.amax(m0 @ y - pyr.fpts[:,1]) < tol
    assert np.amax(m0 @ z - pyr.fpts[:,2]) < tol
    print('Passed M0 test 3: Linear states')


def test_m1():
    pyr = PyrShape(False, cfg)
    m1 = pyr.m1
    tol = 1e-6

    # Test 1: Constant states should give zero grad
    f = np.ones(3*pyr.nupts)
    assert np.amax(m1 @ f) < tol
    print('Passed M1 test 1: Constant')

    # Test 2: Linear x
    x = pyr.upts[:,0]
    y = np.ones_like(x)
    z = np.ones_like(x)
    f = np.array([x,y,z]).reshape(-1, order='C')
    assert np.amax(m1 @ f - 1) < tol

    # Test 2: Linear x
    y = pyr.upts[:,1]
    x = np.ones_like(y)
    z = np.ones_like(y)
    f = np.array([x,y,z]).reshape(-1, order='C')
    assert np.amax(m1 @ f - 1) < tol

    # Test 2: Linear x
    z = pyr.upts[:,2]
    x = np.ones_like(z)
    y = np.ones_like(z)
    f = np.array([x,y,z]).reshape(-1, order='C')
    assert np.amax(m1 @ f - 1) < tol
    print('Passed M1 test 2: Linear')

    if pyr.order > 1:
        # Test 3: Quadratic x
        x = pyr.upts[:,0]
        y = np.ones_like(x)
        z = np.ones_like(x)
        f = np.array([x**2,y,z]).reshape(-1, order='C')
        print(m1 @ f)
        print(2*x)
        assert np.amax(m1 @ f - 2*x) < tol

        # Test 3: Quadratic y
        y = pyr.upts[:,1]
        x = np.ones_like(y)
        z = np.ones_like(y)
        f = np.array([x,y**2,z]).reshape(-1, order='C')
        assert np.amax(m1 @ f - 2*y) < tol

        # Test 2: Quadratic z
        z = pyr.upts[:,2]
        x = np.ones_like(z)
        y = np.ones_like(z)
        f = np.array([x,y,z**2]).reshape(-1, order='C')
        assert np.amax(m1 @ f - 2*z) < tol
        print('Passed M1 test 3: Quadratic')

def test_m3():
    pyr = TetShape(False, cfg)
    m1 = pyr.m1
    m2 = pyr.m2
    m3 = pyr.m3
    tol = 1e-6

    # Test 1: Constant state
    f = np.ones(pyr.nupts*3)
    f2 = np.ones(pyr.nfpts)*np.sum(pyr.norm_fpts, axis=1)

    d = (m1 - m3 @ m2) @ f + m3 @ f2
    assert np.amax(d) < tol
    print('Passed M3 test 1: Constant')

def compute_tri_area_and_norm(x1, x2, x3):
    z1 = x2 - x1
    z2 = x3 - x1

    n = np.cross(z1, z2)
    a = n/2.0
    return n

def test_norms():
    tetl = TetLShape(False, cfg)
    tetr = TetRShape(False, cfg)
    pyr = PyrShape(False, cfg)

    print(pyr.norm_fpts)

    tol = 1e-6

    for ftype, fun, norms in pyr.faces:
        if ftype == 'quad':
            vals = np.array([fun(-1, -1), fun(1, -1), fun(-1, 1), fun(1, 1)])
        elif ftype == 'tri':
            vals = np.array([fun(-1, -1), fun(1, -1), fun(-1, 1)])
        for val in vals:
            assert np.min(val) >= -1-tol and np.max(val) <= 1 + tol
        print(compute_tri_area_and_norm(vals[0], vals[1], vals[2]))

    for ftype, fun, norms in tetl.faces:
        vals = [fun(-1, -1), fun(1, -1), fun(-1, 1)]
        for val in vals:
            assert np.min(val) >= -1-tol and np.max(val) <= 1 + tol

    for ftype, fun, norms in tetr.faces:
        vals = [fun(-1, -1), fun(1, -1), fun(-1, 1)]
        for val in vals:
            assert np.min(val) >= -1-tol and np.max(val) <= 1 + tol
    


''' 
To do:
    Remove linear elements check
    Create quad-fpts using dual tri?
    Create uniform points?


Done:
    Matrices to redo:
    M0: Interpolation to fpts (piece-wise)
    M1: Divergence of interior flux at upts (piece-wise, average at interior diagonal)
    M3: Divergence of interface flux at upts (piece-wise, average at interior diagonal)

    Make map for which upts/fpts indices correspond to which indices of which side half-tet

    Remove interior interface rows from correction matrix?

    Write M0 tests with analytic 
'''

tetl = TetLShape(False, cfg)
tetr = TetRShape(False, cfg)
pyr = PyrShape(False, cfg)

test_m0()
test_m1()
test_m3()
test_norms()