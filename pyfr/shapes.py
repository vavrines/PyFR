# -*- coding: utf-8 -*-

import itertools as it
from functools import cached_property
from math import exp
import re

import numpy as np

from pyfr.nputil import block_diag, clean
from pyfr.polys import get_polybasis
from pyfr.quadrules import get_quadrule


def _proj_pts(projector, pts):
    pts = np.atleast_2d(pts.T)
    return np.vstack(np.broadcast_arrays(*projector(*pts))).T


@clean
def _proj_l2(qrule, basis):
    return basis.vdm.T @ (qrule.wts*basis.ortho_basis_at(qrule.pts))


class BaseShape:
    name = None
    ndims = -1

    npts_coeffs = None
    npts_cdenom = None

    npts_for_face = {
        'line': lambda order: order + 1,
        'tri': lambda order: (order + 1)*(order + 2) // 2,
        'quad': lambda order: (order + 1)**2
    }

    def __init__(self, nspts, cfg):
        self.nspts = nspts
        self.cfg = cfg
        self.order = cfg.getint('solver', 'order')

        self.antialias = cfg.get('solver', 'anti-alias', 'none')
        self.antialias = {s.strip() for s in self.antialias.split(',')}
        self.antialias.discard('none')
        if self.antialias - {'flux', 'surf-flux'}:
            raise ValueError('Invalid anti-alias options')

        self.ubasis = get_polybasis(self.name, self.order + 1, self.upts)

        if nspts:
            self.nsptsord = nsptord = self.order_from_nspts(nspts)
            self.sbasis = get_polybasis(self.name, nsptord, self.spts)

            # Basis for free-stream metric
            # We need p-th order pseudo grid points, which includes
            # p-th order points on faces.
            # It guarantees th q-th order collocation projection on the face
            # on the both adjacent cells.
            # Ref. 1 JCP 281, 28-54, Sec 4.2
            # Ref. 2 JSC 26(3), 301-327, Definition 1
            self.mbasis = get_polybasis(self.name, max(self.order + 1, 2),
                                        self.mpts)

    @classmethod
    def nspts_from_order(cls, sptord):
        return int(np.polyval(cls.npts_coeffs, sptord)) // cls.npts_cdenom

    @classmethod
    def order_from_nspts(cls, nspts):
        # Obtain the coefficients for the poly: P(n) - nspts = 0
        coeffs = list(cls.npts_coeffs)
        coeffs[-1] -= cls.npts_cdenom*int(nspts)

        # Iterate
        for n in range(1, 15):
            if np.polyval(coeffs, n) == 0:
                return n
        else:
            raise ValueError('Invalid number of shape points')

    @clean
    def opmat(self, expr):
        expr = expr.lower().replace('*', '@')

        if not re.match(r'[m0-9\-+@() ]+$', expr):
            raise ValueError('Invalid operator matrix expression')

        mats = {m: getattr(self, m) for m in re.findall(r'm\d+', expr)}
        return eval(expr, {'__builtins__': None}, mats)

    @cached_property
    def m4(self):
        m = self.m1.reshape(self.nupts, -1, self.nupts).swapaxes(0, 1)
        return m.reshape(-1, self.nupts)

    @cached_property
    def m6(self):
        m = self.norm_fpts.T[:, None, :]*self.m3
        return m.reshape(-1, self.nfpts)

    @cached_property
    def m7(self):
        return self.ubasis.nodal_basis_at(self.qpts)

    @cached_property
    def m8(self):
        return _proj_l2(self._eqrule, self.ubasis)

    @property
    def m9(self):
        return block_diag([self.m8]*self.ndims)

    @cached_property
    @clean
    def m10(self):
        ub = self.ubasis

        n = max(sum(dd) for dd in ub.degrees)
        ncut = self.cfg.getint('soln-filter', 'cutoff')
        order = self.cfg.getint('soln-filter', 'order')
        alpha = self.cfg.getfloat('soln-filter', 'alpha')

        A = np.ones(self.nupts)
        for i, d in enumerate(sum(dd) for dd in ub.degrees):
            if d >= ncut < n:
                A[i] = exp(-alpha*((d - ncut)/(n - ncut))**order)

        return np.linalg.solve(ub.vdm, A[:, None]*ub.vdm).T

    @cached_property
    def nupts(self):
        n = self.order + 1
        return int(np.polyval(self.npts_coeffs, n)) // self.npts_cdenom

    @cached_property
    def upts(self):
        if self.name in ['tetl', 'tetr']:
            rname = self.cfg.get(f'solver-elements-tet', 'soln-pts')
            return get_quadrule('tet', rname, self.nupts).pts
        else:
            rname = self.cfg.get(f'solver-elements-{self.name}', 'soln-pts')
            return get_quadrule(self.name, rname, self.nupts).pts
        

    def _get_qrule(self, eleint, kind, **kwargs):
        sect = f'solver-{eleint}-{kind}'

        if self.cfg.hasopt(sect, 'quad-pts'):
            kwargs['rule'] = self.cfg.get(sect, 'quad-pts')

        if self.cfg.hasopt(sect, 'quad-deg'):
            kwargs['qdeg'] = self.cfg.getint(sect, 'quad-deg')

        return get_quadrule(kind, **kwargs)

    @cached_property
    def _eqrule(self):
        return self._get_qrule('elements', self.name)

    @cached_property
    def _iqrules(self):
        return {kind: self._get_qrule('interfaces', kind, flags='s')
                for kind in {k for k, p, n in self.faces}}

    @property
    def qpts(self):
        return self._eqrule.pts

    @property
    def nqpts(self):
        return len(self.qpts)

    @cached_property
    def fpts(self):
        ppts = []

        for kind, proj, norm in self.faces:
            # Obtain the flux points in reference space for the face type
            if 'surf-flux' in self.antialias:
                r = self._iqrules[kind]
            else:
                rule = self.cfg.get(f'solver-interfaces-{kind}', 'flux-pts')
                npts = self.npts_for_face[kind](self.order)

                r = get_quadrule(kind, rule, npts)

            # Project
            ppts.append(_proj_pts(proj, r.pts))

        return np.vstack(ppts)

    @cached_property
    def fpts_wts(self):
        pwts = []

        for kind, proj, norm in self.faces:
            # Obtain the weights in reference space for the face type
            if 'surf-flux' in self.antialias:
                r = self._iqrules[kind]
            else:
                rule = self.cfg.get(f'solver-interfaces-{kind}', 'flux-pts')
                npts = self.npts_for_face[kind](self.order)

                r = get_quadrule(kind, rule, npts)

            pwts.append(r.wts)

        return np.hstack(pwts)

    @cached_property
    def gbasis_coeffs(self):
        coeffs = []

        # Suitable quadrature rules for various face types
        qrule_map = {
            'line': ('gauss-legendre', self.order + 1),
            'quad': ('gauss-legendre', (self.order + 1)**2),
            'tri': ('williams-shunn', 36)
        }

        for kind, proj, norm in self.faces:
            # Obtain a quadrature rule for integrating on the reference face
            # and evaluate this rule at the nodal basis defined by the flux
            # points
            qr = get_quadrule(kind, *qrule_map[kind])
            L = self.facebases[kind].nodal_basis_at(qr.pts)

            # Do the quadrature
            M = self.ubasis.ortho_basis_at(_proj_pts(proj, qr.pts))
            S = np.einsum('i...,ik,ji->kj', qr.wts, L, M)

            coeffs.append(S)

        return np.vstack(coeffs)

    @clean
    def gbasis_at(self, pts):
        return (self.gbasis_coeffs @ self.ubasis.ortho_basis_at(pts)).T

    @property
    def facenorms(self):
        return [norm for kind, proj, norm in self.faces]

    @cached_property
    def norm_fpts(self):
        fnorms = self.facenorms
        return np.vstack([[fn]*n for fn, n in zip(fnorms, self.nfacefpts)])

    @cached_property
    def spts(self):
        return self.std_ele(self.nsptsord - 1)

    @cached_property
    def linspts(self):
        return self.std_ele(1)

    @cached_property
    def facebases(self):
        fb = {}

        for kind in {k for k, p, n in self.faces}:
            rule = self.cfg.get(f'solver-interfaces-{kind}', 'flux-pts')
            npts = self.npts_for_face[kind](self.order)

            pts = get_quadrule(kind, rule, npts).pts

            fb[kind] = get_polybasis(kind, self.order + 1, pts)

        return fb

    @cached_property
    def facefpts(self):
        nf = np.cumsum([0] + self.nfacefpts)
        return [list(range(nf[i], nf[i + 1])) for i in range(len(nf) - 1)]

    @cached_property
    def nfacefpts(self):
        if 'surf-flux' in self.antialias:
            def cnt(k): return len(self._iqrules[k].pts)
        else:
            def cnt(k): return self.npts_for_face[k](self.order)

        return [cnt(kind) for kind, proj, norm in self.faces]

    @property
    def nfpts(self):
        return sum(self.nfacefpts)

    @cached_property
    def mpts(self):
        return self.std_ele(max(self.order, 1))

    @cached_property
    def nmpts(self):
        return len(self.mpts)


class TensorProdShape:
    @classmethod
    def std_ele(cls, sptord):
        pts1d = np.linspace(-1, 1, sptord + 1)
        return [p[::-1] for p in it.product(pts1d, repeat=cls.ndims)]


class QuadShape(TensorProdShape, BaseShape):
    name = 'quad'
    ndims = 2

    # nspts = n^2
    npts_coeffs = [1, 0, 0]
    npts_cdenom = 1

    # Faces: type, reference-to-face projection, normal
    faces = [
        ('line', lambda s: (s, -1), (0, -1)),
        ('line', lambda s: (1, s), (1, 0)),
        ('line', lambda s: (s, 1), (0, 1)),
        ('line', lambda s: (-1, s), (-1, 0)),
    ]

    # Jacobian expressions for a linear element
    jac_exprs = [
        ['((1 - x[1])*V[1][0] - (x[1] + 1)*V[2][0] +'
         ' (x[1] - 1)*V[0][0] + (x[1] + 1)*V[3][0])/4',
         '((1 - x[1])*V[1][1] - (x[1] + 1)*V[2][1] +'
         ' (x[1] - 1)*V[0][1] + (x[1] + 1)*V[3][1])/4'],
        ['((1 - x[0])*V[2][0] - (x[0] + 1)*V[1][0] +'
         ' (x[0] - 1)*V[0][0] + (x[0] + 1)*V[3][0])/4',
         '((1 - x[0])*V[2][1] - (x[0] + 1)*V[1][1] +'
         ' (x[0] - 1)*V[0][1] + (x[0] + 1)*V[3][1])/4']
    ]


class HexShape(TensorProdShape, BaseShape):
    name = 'hex'
    ndims = 3

    # nspts = n^3
    npts_coeffs = [1, 0, 0, 0]
    npts_cdenom = 1

    # Faces: type, reference-to-face projection, normal
    faces = [
        ('quad', lambda s, t: (s, t, -1), (0, 0, -1)),
        ('quad', lambda s, t: (s, -1, t), (0, -1, 0)),
        ('quad', lambda s, t: (1, s, t), (1, 0, 0)),
        ('quad', lambda s, t: (s, 1, t), (0, 1, 0)),
        ('quad', lambda s, t: (-1, s, t), (-1, 0, 0)),
        ('quad', lambda s, t: (s, t, 1), (0, 0, 1)),
    ]

    # Jacobian expressions for a linear element
    jac_exprs = [
        [f'((-x[1]*x[2] + x[1] + x[2] - 1)*V[0][{i}] +'
         f' ( x[1]*x[2] - x[1] - x[2] + 1)*V[1][{i}] +'
         f' ( x[1]*x[2] - x[1] + x[2] - 1)*V[2][{i}] +'
         f' (-x[1]*x[2] + x[1] - x[2] + 1)*V[3][{i}] +'
         f' ( x[1]*x[2] + x[1] - x[2] - 1)*V[4][{i}] +'
         f' (-x[1]*x[2] - x[1] + x[2] + 1)*V[5][{i}] +'
         f' (-x[1]*x[2] - x[1] - x[2] - 1)*V[6][{i}] +'
         f' ( x[1]*x[2] + x[1] + x[2] + 1)*V[7][{i}])/8'
         for i in range(3)],
        [f'((-x[0]*x[2] + x[0] + x[2] - 1)*V[0][{i}] +'
         f' ( x[0]*x[2] - x[0] + x[2] - 1)*V[1][{i}] +'
         f' ( x[0]*x[2] - x[0] - x[2] + 1)*V[2][{i}] +'
         f' (-x[0]*x[2] + x[0] - x[2] + 1)*V[3][{i}] +'
         f' ( x[0]*x[2] + x[0] - x[2] - 1)*V[4][{i}] +'
         f' (-x[0]*x[2] - x[0] - x[2] - 1)*V[5][{i}] +'
         f' (-x[0]*x[2] - x[0] + x[2] + 1)*V[6][{i}] +'
         f' ( x[0]*x[2] + x[0] + x[2] + 1)*V[7][{i}])/8'
         for i in range(3)],
        [f'((-x[0]*x[1] + x[0] + x[1] - 1)*V[0][{i}] +'
         f' ( x[0]*x[1] - x[0] + x[1] - 1)*V[1][{i}] +'
         f' ( x[0]*x[1] + x[0] - x[1] - 1)*V[2][{i}] +'
         f' (-x[0]*x[1] - x[0] - x[1] - 1)*V[3][{i}] +'
         f' ( x[0]*x[1] - x[0] - x[1] + 1)*V[4][{i}] +'
         f' (-x[0]*x[1] + x[0] - x[1] + 1)*V[5][{i}] +'
         f' (-x[0]*x[1] - x[0] + x[1] + 1)*V[6][{i}] +'
         f' ( x[0]*x[1] + x[0] + x[1] + 1)*V[7][{i}])/8'
         for i in range(3)]
    ]
    @cached_property
    def m0(self):
        return self.ubasis.nodal_basis_at(self.fpts)

    @cached_property
    def m1(self):
        m = np.rollaxis(self.ubasis.jac_nodal_basis_at(self.upts), 2)
        return m.reshape(self.nupts, -1)

    @cached_property
    def m2(self):
        m = self.norm_fpts[..., None]*self.m0[:, None, :]
        return m.reshape(self.nfpts, -1)

    @cached_property
    def m3(self):
        m = self.gbasis_at(self.upts)

        if 'surf-flux' in self.antialias:
            fp = [_proj_l2(self._iqrules[kind], self.facebases[kind])
                  for kind, proj, norm in self.faces]

            m = m @ block_diag(fp)

        return m


class TriShape(BaseShape):
    name = 'tri'
    ndims = 2

    # nspts = n*(n + 1)/2
    npts_coeffs = [1, 1, 0]
    npts_cdenom = 2

    # Faces: type, reference-to-face projection, normal
    faces = [
        ('line', lambda s: (s, -1), (0, -1)),
        ('line', lambda s: (-s, s), (1, 1)),
        ('line', lambda s: (-1, s), (-1, 0)),
    ]

    # Jacobian expressions for a linear element
    jac_exprs = [
        [f'(V[{i + 1}][{j}] - V[0][{j}])/2' for j in range(2)]
        for i in range(2)
    ]

    @classmethod
    def std_ele(cls, sptord):
        pts1d = np.linspace(-1, 1, sptord + 1)

        return [(p, q)
                for i, q in enumerate(pts1d)
                for p in pts1d[:(sptord + 1 - i)]]
    @cached_property
    def m0(self):
        return self.ubasis.nodal_basis_at(self.fpts)

    @cached_property
    def m1(self):
        m = np.rollaxis(self.ubasis.jac_nodal_basis_at(self.upts), 2)
        return m.reshape(self.nupts, -1)

    @cached_property
    def m2(self):
        m = self.norm_fpts[..., None]*self.m0[:, None, :]
        return m.reshape(self.nfpts, -1)

    @cached_property
    def m3(self):
        m = self.gbasis_at(self.upts)

        if 'surf-flux' in self.antialias:
            fp = [_proj_l2(self._iqrules[kind], self.facebases[kind])
                  for kind, proj, norm in self.faces]

            m = m @ block_diag(fp)

        return m


class TetShape(BaseShape):
    name = 'tet'
    ndims = 3

    # nspts = n*(n + 1)*(n + 2)/6
    npts_coeffs = [1, 3, 2, 0]
    npts_cdenom = 6

    # Faces: type, reference-to-face projection, normal
    faces = [
        ('tri', lambda s, t: (s, t, -1), (0, 0, -1)),
        ('tri', lambda s, t: (s, -1, t), (0, -1, 0)),
        ('tri', lambda s, t: (-1, t, s), (-1, 0, 0)),
        ('tri', lambda s, t: (s, t, -s - t - 1), (1, 1, 1)),
    ]

    # Jacobian expressions for a linear element
    jac_exprs = [
        [f'(V[{i + 1}][{j}] - V[0][{j}])/2' for j in range(3)]
        for i in range(3)
    ]

    @classmethod
    def std_ele(cls, sptord):
        pts1d = np.linspace(-1, 1, sptord + 1)

        return [(p, q, r)
                for i, r in enumerate(pts1d)
                for j, q in enumerate(pts1d[:(sptord + 1 - i)])
                for p in pts1d[:(sptord + 1 - i - j)]]
    @cached_property
    def m0(self):
        return self.ubasis.nodal_basis_at(self.fpts)

    @cached_property
    def m1(self):
        m = np.rollaxis(self.ubasis.jac_nodal_basis_at(self.upts), 2)
        return m.reshape(self.nupts, -1)

    @cached_property
    def m2(self):
        m = self.norm_fpts[..., None]*self.m0[:, None, :]
        return m.reshape(self.nfpts, -1)

    @cached_property
    def m3(self):
        m = self.gbasis_at(self.upts)

        if 'surf-flux' in self.antialias:
            fp = [_proj_l2(self._iqrules[kind], self.facebases[kind])
                  for kind, proj, norm in self.faces]

            m = m @ block_diag(fp)

        return m


class PriShape(BaseShape):
    name = 'pri'
    ndims = 3

    # nspts = n^2*(n + 1)/2
    npts_coeffs = [1, 1, 0, 0]
    npts_cdenom = 2

    # Faces: type, reference-to-face projection, normal
    faces = [
        ('tri', lambda s, t: (s, t, -1), (0, 0, -1)),
        ('tri', lambda s, t: (s, t, 1), (0, 0, 1)),
        ('quad', lambda s, t: (s, -1, t), (0, -1, 0)),
        ('quad', lambda s, t: (-s, s, t), (1, 1, 0)),
        ('quad', lambda s, t: (-1, s, t), (-1, 0, 0)),
    ]

    # Jacobian expressions for a linear element
    _jac_exprs_xy = [
        [f'((x[2] - 1)*V[0][{j}] + (1 - x[2])*V[{i + 1}][{j}] -'
         f' (x[2] + 1)*V[3][{j}] + (x[2] + 1)*V[{i + 4}][{j}])/4'
         for j in range(3)]
        for i in range(2)
    ]
    _jac_exprs_z = [[f'((x[0] + 1)*V[4][{j}] + (x[0] + x[1])*V[0][{j}] -'
                     f' (x[0] + 1)*V[1][{j}] - (x[0] + x[1])*V[3][{j}] -'
                     f' (x[1] + 1)*V[2][{j}] + (x[1] +    1)*V[5][{j}])/4'
                     for j in range(3)]]
    jac_exprs = _jac_exprs_xy + _jac_exprs_z


    @classmethod
    def std_ele(cls, sptord):
        pts1d = np.linspace(-1, 1, sptord + 1)

        return [(p, q, r)
                for r in pts1d
                for i, q in enumerate(pts1d)
                for p in pts1d[:(sptord + 1 - i)]]


class TetLShape(BaseShape):
    '''
    This tet defined by points
    
    (-1, -1, -1)
    ( 1, -1, -1)
    (-1,  1, -1)
    ( 0,  0,  1)

    Faces:
    bottom
    south
    west
    north-east

    ''' 

    name = 'tetl'
    ndims = 3

    # nspts = n*(n + 1)*(n + 2)/6
    npts_coeffs = [1, 3, 2, 0]
    npts_cdenom = 6

    # Faces: type, reference-to-face projection, normal

    # CHANGE THIS TO NEW TET


    faces = [
        ('tri', lambda s, t: (s, t, -1), (0, 0, -1)),
        ('tri', lambda s, t: (s + (t + 1)/2, (t - 1)/2, t), (0, -1, 0.5)),
        ('tri', lambda s, t: ((s - 1)/2, t + (s + 1)/2, s), (-1, 0, 0.5)),
        ('tri', lambda s, t: (s - (s+t)/2, t - (s+t)/2, -(s+t) - 1), (1, 1, 0))
    ]


    def __init__(self, nspts, cfg):
        super().__init__(nspts, cfg)
        # Move upper sol point to (0,0,1) (instead of (-1, -1, 1))
        rfac = 0.5*(self.upts[:,2] + 1.)
        self.upts[:,0] += rfac
        self.upts[:,1] += rfac
        self.ubasis = get_polybasis('tetl', self.order + 1, self.upts)

    @classmethod
    def std_ele(cls, sptord):
        pts1d = np.linspace(-1, 1, sptord + 1)

        vals = [(p, q, r)
                for i, r in enumerate(pts1d)
                for j, q in enumerate(pts1d[:(sptord + 1 - i)])
                for p in pts1d[:(sptord + 1 - i - j)]]

        # Transform to new tet (top point moves from (-1,-1) tp (0,0))
        newvals = []
        for (p,q,r) in vals:
            rfac = (r + 1.0)/2.0
            pnew = p + rfac
            qnew = q + rfac
            newvals.append( (pnew, qnew, r))

        return newvals
    @cached_property
    def m0(self):
        return self.ubasis.nodal_basis_at(self.fpts)

    @cached_property
    def m1(self):
        m = np.rollaxis(self.ubasis.jac_nodal_basis_at(self.upts), 2)
        return m.reshape(self.nupts, -1)

    @cached_property
    def m2(self):
        m = self.norm_fpts[..., None]*self.m0[:, None, :]
        return m.reshape(self.nfpts, -1)

    @cached_property
    def m3(self):
        m = self.gbasis_at(self.upts)

        if 'surf-flux' in self.antialias:
            fp = [_proj_l2(self._iqrules[kind], self.facebases[kind])
                  for kind, proj, norm in self.faces]

            m = m @ block_diag(fp)

        return m

class TetRShape(BaseShape):

    '''
    This tet defined by points (TetL rotated 180 deg)
    
    ( 1,  1)
    (-1,  1)
    ( 1, -1)
    ( 0,  0)

    Faces:
    bottom
    north
    east
    south-west
    ''' 

    name = 'tetr'
    ndims = 3

    # nspts = n*(n + 1)*(n + 2)/6
    npts_coeffs = [1, 3, 2, 0]
    npts_cdenom = 6

    # Faces: type, reference-to-face projection, normal

    faces = [
        ('tri', lambda s, t: (-s, -t, -1), (0, 0, -1)),
        ('tri', lambda s, t: (-(s + (t + 1)/2), -((t - 1)/2), t), (0, -1, 0.5)),
        ('tri', lambda s, t: (-((s - 1)/2), -(t + (s + 1)/2), s), (-1, 0, 0.5)),
        ('tri', lambda s, t: (-(s - (s+t)/2), -(t - (s+t)/2), -(s+t) - 1), (1, 1, 0))
    ]

    def __init__(self, nspts, cfg):
        super().__init__(nspts, cfg)
        # Move upper sol point to (0,0,1) (instead of (-1, -1, 1))
        rfac = 0.5*(self.upts[:,2] + 1.)
        self.upts[:,0] += rfac
        self.upts[:,1] += rfac
        # Flip upts/fpts
        self.upts[:,:2] *= -1
        # self.fpts[:,:2] *= -1
        self.ubasis = get_polybasis('tetr', self.order + 1, self.upts)

    @classmethod
    def std_ele(cls, sptord):
        pts1d = np.linspace(-1, 1, sptord + 1)

        vals = [(p, q, r)
                for i, r in enumerate(pts1d)
                for j, q in enumerate(pts1d[:(sptord + 1 - i)])
                for p in pts1d[:(sptord + 1 - i - j)]]

        # Transform to new tet (top point moves from (-1,-1) tp (0,0))
        # and rotate 180 deg
        newvals = []
        for (p,q,r) in vals:
            rfac = (r + 1.0)/2.0
            pnew = p + rfac
            qnew = q + rfac
            newvals.append( (-pnew, -qnew, r)) # rotate

        return newvals
    @cached_property
    def m0(self):
        return self.ubasis.nodal_basis_at(self.fpts)

    @cached_property
    def m1(self):
        m = np.rollaxis(self.ubasis.jac_nodal_basis_at(self.upts), 2)
        return m.reshape(self.nupts, -1)

    @cached_property
    def m2(self):
        m = self.norm_fpts[..., None]*self.m0[:, None, :]
        return m.reshape(self.nfpts, -1)

    @cached_property
    def m3(self):
        m = self.gbasis_at(self.upts)

        if 'surf-flux' in self.antialias:
            fp = [_proj_l2(self._iqrules[kind], self.facebases[kind])
                  for kind, proj, norm in self.faces]

            m = m @ block_diag(fp)

        return m

class PyrShape(BaseShape):
    name = 'pyr'
    ndims = 3



    # nspts = n*(n + 1)*(2*n + 1)/6
    npts_coeffs = [2, 3, 1, 0]
    npts_cdenom = 6

    # Faces: type, reference-to-face projection, normal
    '''
    Faces:
    0: bottom
    1: south
    2: east
    3: north
    4: west

    '''
    faces = [
        ('quad', lambda s, t: (s, t, -1), (0, 0, -1)),
        ('tri', lambda s, t: (s + (t + 1)/2, (t - 1)/2, t), (0, -1, 0.5)),
        ('tri', lambda s, t: ((1 - t)/2, -s - (t + 1)/2, t), (1, 0, 0.5)),
        ('tri', lambda s, t: (-s - (t + 1)/2, (1 - t)/2, t), (0, 1, 0.5)),
        ('tri', lambda s, t: ((t - 1)/2, s + (t + 1)/2, t), (-1, 0, 0.5)),
    ]

    # Jacobian expressions for a linear element
    jac_exprs = [
        ['((1 - x[1])*V[1][0] + (-x[1] - 1)*V[2][0] +'
         ' (x[1] - 1)*V[0][0] + (x[1] + 1)*V[3][0])/4',
         '((1 - x[1])*V[1][1] + (-x[1] - 1)*V[2][1] +'
         '(x[1] - 1)*V[0][1] + (x[1] + 1)*V[3][1])/4',
         '((1 - x[1])*V[1][2] + (-x[1] - 1)*V[2][2] +'
         ' (x[1] - 1)*V[0][2] + (x[1] + 1)*V[3][2])/4'],
        ['((1 - x[0])*V[2][0] + (-x[0] - 1)*V[1][0] +'
         ' (x[0] - 1)*V[0][0] + (x[0] + 1)*V[3][0])/4',
         '((1 - x[0])*V[2][1] + (-x[0] - 1)*V[1][1] +'
         ' (x[0] - 1)*V[0][1] + (x[0] + 1)*V[3][1])/4',
         '((1 - x[0])*V[2][2] + (-x[0] - 1)*V[1][2] +'
         '(x[0] - 1)*V[0][2] + (x[0] + 1)*V[3][2])/4'],
        ['(-V[0][0] - V[1][0] - V[2][0] - V[3][0] + 4*V[4][0])/8',
         '(-V[0][1] - V[1][1] - V[2][1] - V[3][1] + 4*V[4][1])/8',
         '(-V[0][2] - V[1][2] - V[2][2] - V[3][2] + 4*V[4][2])/8']
    ]



    def __init__(self, nspts, cfg):
        super().__init__(nspts, cfg)
        self.tetl = TetLShape(False, self.cfg)
        self.tetr = TetRShape(False, self.cfg)

        p = self.order
        self.ntripts = (p+1)*(p+2)//2
        self.nquadpts = (p+1)**2

        self.make_fpts_map()
        self.make_upts_map()


    @classmethod
    def std_ele(cls, sptord):
        npts1d = 2*sptord + 1
        pts1d = np.linspace(-1, 1, npts1d)

        return [(p, q, r)
                for i, r in enumerate(pts1d[::2])
                for q in pts1d[i:npts1d - i:2]
                for p in pts1d[i:npts1d - i:2]]

    def make_fpts_map(self):
        self.fpts_map = [None]*self.nfpts
        self.fpts_idxs = [None]*self.nfpts
        n = 0
        # Bottom quad points
        # Format: i/x changes first, then j/y

        tetl_bot = self.tetl.fpts[:self.ntripts]
        tetr_bot = self.tetr.fpts[:self.ntripts]

        tol = 1e-8
        for i in range(self.nquadpts):
            coord = self.fpts[n, :]
            [x,y,z] = coord

            if np.abs(x+y) < tol:
                tag = 'm'
                idxl = np.where((np.abs(tetl_bot - coord) < tol).all(axis=1))[0][0]
                idxr = np.where((np.abs(tetr_bot - coord) < tol).all(axis=1))[0][0]
                # assert idxl == idxr, f'{idxl}, {idxr}'
                self.fpts_idxs[n] = (idxl, idxr)
            elif x + y < -tol:
                tag = 'l'
                self.fpts_idxs[n] = np.where((np.abs(tetl_bot - coord) < tol).all(axis=1))[0][0]
            else:
                tag = 'r'
                self.fpts_idxs[n] = np.where((np.abs(tetr_bot - coord) < tol).all(axis=1))[0][0]

            self.fpts_map[n] = tag

            n += 1


        '''
        
        Pyr faces:
            0: bottom
            1: south
            2: east
            3: north
            4: west

        TetL faces:
            0: bottom
            1: south
            2: west
            3: north-east

        TetR faces:
            0: bottom
            1: north
            2: east
            3: south-west
        
        '''
        for face in range(1,5):
            pyr_coords = self.fpts[self.nquadpts + (face-1)*self.ntripts: \
                                   self.nquadpts + (face)*self.ntripts, : ]
            if face == 1:
                tag = 'l'
                offset = 1*self.ntripts
                tet_coords = self.tetl.fpts[offset:offset+self.ntripts,:]
            elif face == 2:
                tag = 'r'
                offset = 2*self.ntripts
                tet_coords = self.tetr.fpts[offset:offset+self.ntripts,:]
            elif face == 3:
                tag = 'r'
                offset = 1*self.ntripts
                tet_coords = self.tetr.fpts[offset:offset+self.ntripts,:]
            elif face == 4:
                tag = 'l'
                offset = 2*self.ntripts
                tet_coords = self.tetl.fpts[offset:offset+self.ntripts,:]

            for i in range(self.ntripts):
                self.fpts_map[n] = tag
                coord = self.fpts[n, :]

                locidx = np.where((tet_coords == coord).all(axis=1))[0][0]
                self.fpts_idxs[n] = offset + locidx
                n += 1

    def make_upts_map(self):
        self.upts_map = [None]*self.nupts
        self.upts_idxs = [None]*self.nupts

        tol = 1e-8
        for i in range(self.nupts):
            coord = self.upts[i, :]
            (x,y,z) = coord

            if np.abs(x+y) < tol: # Middle/diag points
                lidx = np.where((self.tetl.upts == coord).all(axis=1))[0][0]
                ridx = np.where((self.tetr.upts == coord).all(axis=1))[0][0]
                self.upts_map[i] = 'm'
                self.upts_idxs[i] = (lidx, ridx)
            elif x+y > tol:
                self.upts_map[i] = 'r'
                idx = np.where((self.tetr.upts == coord).all(axis=1))[0][0]
                self.upts_idxs[i] = idx
            elif x+y < -tol:
                self.upts_map[i] = 'l'
                idx = np.where((self.tetl.upts == coord).all(axis=1))[0][0]
                self.upts_idxs[i] = idx

    def expand_upts_row(self, row, side):
        newrow = np.zeros(self.nupts)
        #Arbitrarily choose right side to include midpoints as well
        for i,v in enumerate(row):
            newidx = None
            if side == 'right':
                for j, vv in enumerate(self.upts_idxs):
                    if isinstance(vv, np.int64):
                        if vv == i and self.upts_map[j] == 'r':
                            newidx = j
                    elif isinstance(vv, tuple):
                        if vv[1] == i:
                            newidx = j
            elif side == 'left':
                for j, vv in enumerate(self.upts_idxs):
                    if isinstance(vv, np.int64):
                        if vv == i and self.upts_map[j] == 'l':
                            newidx = j
                    elif isinstance(vv, tuple):
                        if vv[0] == i:
                            newidx = j

            assert newidx is not None
            newrow[newidx] = v
        return newrow

    def expand_fpts_row(self, row, side):
        newrow = np.ones(self.nfpts)*np.nan
        #Arbitrarily choose right side to include midpoints as well
        for i,v in enumerate(row):
            newidx = None
            if side == 'right':
                for j, vv in enumerate(self.fpts_idxs):
                    if isinstance(vv, np.int64):
                        if vv == i and self.fpts_map[j] == 'r':
                            newidx = j
                    elif isinstance(vv, tuple):
                        if vv[1] == i:
                            newidx = j
            elif side == 'left':
                for j, vv in enumerate(self.fpts_idxs):
                    if isinstance(vv, np.int64):
                        if vv == i and self.fpts_map[j] == 'l':
                            newidx = j
                    elif isinstance(vv, tuple):
                        if vv[0] == i:
                            newidx = j

            assert newidx is not None
            newrow[newidx] = v
        return newrow
    
    @cached_property
    def m0(self):
        M = np.zeros((self.nfpts, self.nupts))
        tol = 1e-8

        for i, (x,y,z) in enumerate(self.fpts):
            if np.abs(x + y) < tol: # Midpoints
                rowl = self.tetl.ubasis.nodal_basis_at([[x,y,z]])[0]
                rowr = self.tetr.ubasis.nodal_basis_at([[x,y,z]])[0]
                M[i, :] =  0.5*self.expand_upts_row(rowl, 'left')
                M[i, :] += 0.5*self.expand_upts_row(rowr, 'right')
            elif (x + y) > tol: # Right
                row = self.tetr.ubasis.nodal_basis_at([[x,y,z]])[0]
                M[i, :] = self.expand_upts_row(row, 'right')
            else: # Left
                row = self.tetl.ubasis.nodal_basis_at([[x,y,z]])[0]
                M[i, :] = self.expand_upts_row(row, 'left')

        return M

    @cached_property
    def m1(self):
        M = np.zeros((self.nupts, 3*self.nupts))
        tol = 1e-8

        def expand_grad(mat, side):			
            r1 = self.expand_upts_row(mat[0,0,:], side)
            r2 = self.expand_upts_row(mat[0,1,:], side)
            r3 = self.expand_upts_row(mat[0,2,:], side)

            newmat = np.array([r1, r2, r3])
            return newmat

        for i, (x,y,z) in enumerate(self.upts):
            if np.abs(x + y) < tol: # Midpoints
                ml = np.rollaxis(self.tetl.ubasis.jac_nodal_basis_at([[x,y,z]]), 2)
                ml = expand_grad(ml, 'left')

                mr = np.rollaxis(self.tetr.ubasis.jac_nodal_basis_at([[x,y,z]]), 2)
                mr = expand_grad(mr, 'right')

                M[i, :] =  0.5*ml.reshape(1, -1) + 0.5*mr.reshape(1, -1)
            elif (x + y) > tol: # Right
                m = np.rollaxis(self.tetr.ubasis.jac_nodal_basis_at([[x,y,z]]), 2)
                m = expand_grad(m, 'right')
                M[i, :] = m.reshape(1, -1)
            else: # Left 
                m = np.rollaxis(self.tetl.ubasis.jac_nodal_basis_at([[x,y,z]]), 2)
                m = expand_grad(m, 'left')
                M[i, :] = m.reshape(1, -1)

        return M

    @cached_property
    def m2(self):
        m = self.norm_fpts[..., None]*self.m0[:, None, :]
        return m.reshape(self.nfpts, -1)

    @cached_property
    def m3(self):
        M = np.zeros((self.nupts, self.nfpts))
        tol = 1e-8

        # Remove contribution of interior interface (tet diagonal)
        def remove_diag(m):
            return m[:-self.ntripts]
        
        def average_rows(r1, r2):
            r3 = np.zeros_like(r1)
            for i in range(len(r1)):
                if np.isnan(r1[i]):
                    r3[i] = r2[i]
                elif np.isnan(r2[i]):
                    r3[i] = r1[i]
                else:
                    r3[i] = 0.5*(r1[i] + r2[i])
            return r3
        
                
        for i, (x,y,z) in enumerate(self.upts):
            ml = remove_diag(self.tetl.gbasis_at([[x,y,z]])[0])
            mr = remove_diag(self.tetr.gbasis_at([[x,y,z]])[0])

            if np.abs(x + y) < tol: # Midpoints
                M[i, :] = average_rows(self.expand_fpts_row(ml, 'left'), self.expand_fpts_row(mr, 'right'))
            elif (x + y) > tol: # Right
                M[i, :] = np.nan_to_num(self.expand_fpts_row(mr, 'right'), 0.0)
            else: # Left 
                M[i, :] = np.nan_to_num(self.expand_fpts_row(ml, 'left'), 0.0)

        return M
