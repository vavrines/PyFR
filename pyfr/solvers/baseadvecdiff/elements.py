# -*- coding: utf-8 -*-

from pyfr.backends.base.kernels import ComputeMetaKernel
from pyfr.solvers.baseadvec import BaseAdvectionElements
from pyfr.quadrules import get_quadrule
import numpy as np

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
			lambdatype = self.cfg.get('solver', 'lambda_type', 'davis')

			# Template arguments
			tplargs = dict(
				nvars=self.nvars, nupts=self.nupts, ndims=self.ndims, svar=shockvar,
				cav=self.cfg.items_as('solver-artificial-viscosity', float),
				c=self.cfg.items_as('constants', float),
				order=self.basis.order, ubdegs=ubdegs,
				invvdm=self.basis.ubasis.invvdm.T,
				adjmat=adjmat, lambda_type=lambdatype,
				weights=weights
			)



			# Allocate space for the artificial viscosity vector
			self.artvisc = backend.matrix((self.nupts, self.nvars, self.neles),
										  extent=nonce + 'artvisc', tags=tags)


			# Apply the sensor to estimate the required artificial viscosity
			self.kernels['shocksensor'] = lambda: backend.kernel(
				'shocksensor', tplargs=tplargs, dims=[self.neles],
				u=self.scal_upts_inb, uf=self._scal_fpts, 
				plocu=self.ploc_at('upts'), plocf=self.ploc_at('fpts'),
				artvisc=self.artvisc, rcpdjac=self.rcpdjac_at('upts')
			)

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
		rule = self.cfg.get('solver-elements-hex', 'soln-pts', 'gauss-legendre')

		rule = get_quadrule(eletype, rule=None, npts=nupts, qdeg=p, flags=None)
		return rule.wts