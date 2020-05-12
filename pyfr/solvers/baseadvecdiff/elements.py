# -*- coding: utf-8 -*-

from pyfr.backends.base.kernels import ComputeMetaKernel
from pyfr.solvers.baseadvec import BaseAdvectionElements
from pyfr.quadrules import get_quadrule
import numpy as np
from scipy import interpolate

class BaseAdvectionDiffusionElements(BaseAdvectionElements):
	@property
	def _scratch_bufs(self):
		bufs = {'scal_fpts', 'vect_fpts', 'vect_upts'}

		if 'flux' in self.antialias:
			bufs |= {'scal_qpts', 'vect_qpts'}
		elif 'div-flux' in self.antialias:
			bufs |= {'scal_qpts'}

		if self._soln_in_src_exprs:
			if 'div-flux' in self.antialias:
				bufs |= {'scal_qpts_cpy'}
			else:
				bufs |= {'scal_upts_cpy'}

		bufs |= {'scal_fpts_cpy'}

		return bufs

	def set_backend(self, backend, nscalupts, nonce):
		super().set_backend(backend, nscalupts, nonce)

		# Register pointwise kernels
		backend.pointwise.register(
			'pyfr.solvers.baseadvecdiff.kernels.gradcoru'
		)

		self.kernels['_copy_fpts'] = lambda: backend.kernel(
			'copy', self._vect_fpts.rslice(0, self.nfpts), self._scal_fpts
		)
		self.kernels['tgradpcoru_upts'] = lambda: backend.kernel(
			'mul', self.opmat('M4 - M6*M0'), self.scal_upts_inb,
			out=self._vect_upts
		)
		self.kernels['tgradcoru_upts'] = lambda: backend.kernel(
			'mul', self.opmat('M6'), self._vect_fpts.rslice(0, self.nfpts),
			 out=self._vect_upts, beta=1.0
		)
		self.kernels['gradcoru_upts'] = lambda: backend.kernel(
			'gradcoru', tplargs=dict(ndims=self.ndims, nvars=self.nvars),
			 dims=[self.nupts, self.neles], smats=self.smat_at('upts'),
			 rcpdjac=self.rcpdjac_at('upts'), gradu=self._vect_upts
		)

		def gradcoru_fpts():
			nupts, nfpts = self.nupts, self.nfpts
			vupts, vfpts = self._vect_upts, self._vect_fpts

			# Exploit the block-diagonal form of the operator
			muls = [backend.kernel('mul', self.opmat('M0'),
								   vupts.rslice(i*nupts, (i + 1)*nupts),
								   vfpts.rslice(i*nfpts, (i + 1)*nfpts))
					for i in range(self.ndims)]

			return ComputeMetaKernel(muls)

		self.kernels['gradcoru_fpts'] = gradcoru_fpts

		if 'flux' in self.antialias:
			def gradcoru_qpts():
				nupts, nqpts = self.nupts, self.nqpts
				vupts, vqpts = self._vect_upts, self._vect_qpts

				# Exploit the block-diagonal form of the operator
				muls = [backend.kernel('mul', self.opmat('M7'),
									   vupts.rslice(i*nupts, (i + 1)*nupts),
									   vqpts.rslice(i*nqpts, (i + 1)*nqpts))
						for i in range(self.ndims)]

				return ComputeMetaKernel(muls)

			self.kernels['gradcoru_qpts'] = gradcoru_qpts

		# Shock capturing
		shock_capturing = self.cfg.get('solver', 'shock-capturing', 'none')
		if shock_capturing == 'artificial-viscosity':
			tags = {'align'}

			# Register the kernels
			backend.pointwise.register(
				'pyfr.solvers.baseadvecdiff.kernels.shocksensor'
			)

			# Obtain the scalar variable to be used for shock sensing
			shockvar = self.convarmap[self.ndims].index(self.shockvar)

			# Obtain the degrees of the polynomial modes in the basis
			ubdegs = [sum(dd) for dd in self.basis.ubasis.degrees]

			adjmat = self.adj_mat(self.basis.upts, self.basis.fpts)
			weights = self.get_weights(self.basis.upts)
			lambdatype = self.cfg.get('artificial-viscosity', 'lambda_type', 'davis')
			gradmat = self.grad_matrices(self.basis.upts, self.basis.fpts, adjmat)
			ffaces = self.flux_faces() # Face index corresponding to flux point
			prefac = self.calc_prefactor()


			# Template arguments
			tplargs = dict(
				nvars=self.nvars, nupts=self.nupts, nfpts=self.nfpts, ndims=self.ndims,
				cav=self.cfg.items_as('solver-artificial-viscosity', float),
				c=self.cfg.items_as('constants', float),
				order=self.basis.order, ubdegs=ubdegs,
				invvdm=self.basis.ubasis.invvdm.T,
				adjmat=adjmat, lambda_type=lambdatype,
				weights=weights, gradmat=gradmat, ffaces=ffaces,
				prefac=prefac
			)



			# Allocate space for the artificial viscosity vector
			self.artvisc = backend.matrix((self.nupts, self.nvars, self.neles),
										  extent=nonce + 'artvisc', tags=tags)


			# Apply the sensor to estimate the required artificial viscosity
			self.kernels['shocksensor'] = lambda: backend.kernel(
				'shocksensor', tplargs=tplargs, dims=[self.neles],
				u=self.scal_upts_inb, uf=self._scal_fpts, 
				artvisc=self.artvisc, urcpdjac=self.rcpdjac_at('upts'), frcpdjac=self.rcpdjac_at('fpts'),
				usmats=self.ele_smat_at('upts'), rcpfsmats=self.rcp_ele_smat_at('fpts')
			)

			# ELE_SMATS FORMAT:
	        # [0] = dxi/dx,   [1] = dxi/dy,   [2] = dxi/dz
	        # [3] = deta/dx,  [4] = deta/dy,  [5] = deta/dz
	        # [6] = dzeta/dx, [7] = dzeta/dy, [8] = dzeta/dz

			# Add as source term in negdivconf
			backend.pointwise.register('pyfr.solvers.baseadvecdiff.kernels.negdivconfav')
			plocupts = self.ploc_at('upts') if self._ploc_in_src_exprs else None
			solnupts = self._scal_upts_cpy if self._soln_in_src_exprs else None 
				   
			srctplargs = {
				'ndims': self.ndims,
				'nvars': self.nvars,
				'srcex': self._src_exprs
			}

			if self._ploc_in_src_exprs:
				self.kernels['copy_soln'] = lambda: backend.kernel(
					'copy', self._scal_upts_cpy, self.scal_upts_inb
				)

			self.kernels['negdivconf'] = lambda: backend.kernel(
				'negdivconfav', tplargs=srctplargs,
				dims=[self.nupts, self.neles], tdivtconf=self.scal_upts_outb,
				rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=solnupts,
				artvisc=self.artvisc
			)
		elif shock_capturing == 'none':
			self.artvisc = None
		else:
			raise ValueError('Invalid shock capturing scheme')

	def get_artvisc_fpts_for_inter(self, eidx, fidx):        
		nfp = self.nfacefpts[fidx]
		rmap = self._srtd_face_fpts[fidx][eidx]
		cmap = (eidx,)*nfp

		return (self.artvisc.mid,)*nfp, rmap, cmap

	def adj_mat(self, upts, fpts):
		p = self.basis.order
		nupts = len(upts)

		uidx = lambda xidx, yidx, zidx: xidx + yidx*(p+1) + zidx*((p+1)**2)
		cartidx = lambda pidx: (pidx%(p+1), (pidx%(p+1)**2)//(p+1), pidx//((p+1)**2))

		def findidx(xidx, yidx, zidx):
			face = None
			if xidx == -1:
				face = 4
			elif xidx == p+1:
				face = 2
			if yidx == -1:
				face = 1
			elif yidx == p+1:
				face = 3
			if zidx == -1:
				face = 0
			elif zidx == p+1:
				face = 5

			if face == None:
				return uidx(xidx, yidx, zidx)

			if face == 4 or face == 2:
				return (face*(p+1)**2 + yidx + zidx*(p+1)) + nupts

			if face == 1 or face == 3:
				return (face*(p+1)**2 + xidx + zidx*(p+1)) + nupts

			if face == 0 or face == 5:
				return (face*(p+1)**2 + xidx + yidx*(p+1)) + nupts


		M = np.zeros((nupts, 6))
		# [-x, +x, -y, +y, -z, +z] neighbors
		for pidx in range(nupts):
			(xidx, yidx, zidx) = cartidx(pidx)
			M[pidx,0] = findidx(xidx-1, yidx, zidx)
			M[pidx,1] = findidx(xidx+1, yidx, zidx)
			M[pidx,2] = findidx(xidx, yidx-1, zidx)
			M[pidx,3] = findidx(xidx, yidx+1, zidx)
			M[pidx,4] = findidx(xidx, yidx, zidx-1)
			M[pidx,5] = findidx(xidx, yidx, zidx+1)

		return M.astype('i4')

	def get_weights(self, upts):
		p = self.basis.order
		nupts = len(upts)
		eletype = 'hex'
		#rule = self.cfg.get('solver-elements-hex', 'soln-pts', 'gauss-legendre')
		rule = 'gauss-legendre'

		qrule = get_quadrule(eletype, rule=rule, npts=nupts, qdeg=p, flags=None)
		return qrule.wts



	def grad_matrices(self, upts3d, fpts3d, adj_mat):		

		p = self.basis.order
		upts1d = get_quadrule('line', rule='gauss-legendre', npts=(p+1), qdeg=p, flags=None).pts
		nupts = len(upts3d)
		nfpts = len(fpts3d)

		GM = np.zeros((nupts, 6, 3))

		def Lag1d(pts, idx):
			vals = np.zeros(len(pts))
			vals[idx] = 1.
			return interpolate.lagrange(pts, vals)

		derivLags1d = []
		for idx in range(len(upts1d)):
			lag1d = Lag1d(upts1d, idx)
			derivLags1d.append(lag1d.deriv())

		# Returns the gradient of the j-th (jx, jy, jz) interpolating polynomial (in 3D space) at point (x,y,z)
		def gradLag(upts1d, jx, jy, jz, x,y,z):
			dLx = derivLags1d[jx]
			dLy = derivLags1d[jy]
			dLz = derivLags1d[jz]
			grad = [dLx(x), dLy(y), dLz(z)]
			return np.array(grad)

		# Returns the gradient of the correction function at the interface (in 1D)
		def gradCorrFunc(p, deltaxi):
			c = np.zeros(p+2)
			c[p] = 0.5
			c[p+1] = 0.5
			L = np.polynomial.legendre.Legendre(c)
			Ld = L.deriv()
			if deltaxi < 0. or deltaxi > 1.:
				raise ValueError('Deltaxi out of bounds [0,1]:', dletaxi)
			return Ld(1.0 - deltaxi)




		cartidx = lambda pidx: (pidx%(p+1), (pidx%(p+1)**2)//(p+1), pidx//((p+1)**2)) # Global solution idx (pidx) to xidx, yidx, zidx
		face    = lambda fidx: fidx//((p+1)**2) # Face number from fpts idx

		# Calculate distance from flux point to solution point in computational space (kind of a hack)
		deltaxi = 1.0 - max(upts1d)

		# HEX FACE ORDERING:
		# 0 -> z = -1
		# 1 -> y = -1
		# 2 -> x =  1
		# 3 -> y =  1
		# 4 -> x = -1
		# 5 -> z =  1
		def gradFluxPoint(fidx, deltaxi):		
			dcorr = gradCorrFunc(p, deltaxi)
			grads = [[0,     0,     -dcorr],  # z = -1
					 [0,     -dcorr,0     ],  # y = -1
					 [dcorr, 0,     0     ],  # x =  1
					 [0,     dcorr, 0     ],  # y =  1
					 [-dcorr,0,     0     ],  # x = -1
					 [0,     0,     dcorr ],] # z =  1

			return grads[face(fidx)]


		for pidx in range(nupts):
			for j in range(6):
				adjidx = adj_mat[pidx][j]

				# If adjacent point is a solution point
				if adjidx < nupts:
					# Get tensor product position indices
					(jx, jy, jz) = cartidx(adjidx)			
					# Evaluate at solution point pidx
					[x,y,z] = upts3d[pidx]
					grad = gradLag(upts1d, jx, jy, jz, x,y,z)

				else:
					fidx = adjidx - nupts
					grad = gradFluxPoint(fidx, deltaxi)
				GM[pidx, j, :] = grad	
		return GM

	def flux_faces(self):	# Face index corresponding to flux point
		# HEX FACE ORDERING:
		# 0 -> z = -1
		# 1 -> y = -1
		# 2 -> x =  1
		# 3 -> y =  1
		# 4 -> x = -1
		# 5 -> z =  1
		p = self.basis.order
		nfpts = self.nfpts
		ffaces = np.zeros(nfpts, dtype=int)
		for fidx in range(nfpts):
			ffaces[fidx] = int(fidx//((p+1)**2))
		return ffaces

	def calc_prefactor(self): 
		# prefactor = 0.5*phi_i(x_j)*w_j'/w_i evaluated on flux point x_j for solution basis function phi_i if using only adjacent points (so Phi_i is first/last GL point)
		# Weights correspond to 2d (w') and 3d (w) quad weights
		p = self.basis.order
		upts1d = get_quadrule('line', rule='gauss-legendre', npts=(p+1), qdeg=p, flags=None).pts
		w1d = get_quadrule('line', rule='gauss-legendre', npts=(p+1), qdeg=p, flags=None).wts		
		vals = np.zeros(len(upts1d))
		vals[0] = 1.
		L = interpolate.lagrange(upts1d, vals)

		phi_val = L(-1.0)
		w = w1d[0] # Since w_j' = 2D weight and w_i is 3D weight along the same 2d indices, just one point off in 3d
		return 0.5*phi_val/w


