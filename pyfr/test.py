from pyfr.shapes import *
from pyfr.polys import *
from pyfr.inifile import Inifile



path = '/mnt/c/Users/Tarik/Documents/GitHub/PyFR-PANS/pyfr/test.ini'
cfg = Inifile.load(path)


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

	name = 'tet'
	ndims = 3

	# nspts = n*(n + 1)*(n + 2)/6
	npts_coeffs = [1, 3, 2, 0]
	npts_cdenom = 6

	# Faces: type, reference-to-face projection, normal

	# CHANGE THIS TO NEW TET


	faces = [
		('tri', lambda s, t: (s, t, -1), (0, 0, -1)),
		('tri', lambda s, t: (s + (t + 1)/2, (t - 1)/2, t), (0, -1, 0.5)),
		('tri', lambda s, t: ((t - 1)/2, s + (t + 1)/2, t), (-1, 0, 0.5)),
		('tri', lambda s, t: (-s - (t + 1)/2, s + (t + 1)/2, t), (1, 1, 0))
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

	name = 'tet'
	ndims = 3

	# nspts = n*(n + 1)*(n + 2)/6
	npts_coeffs = [1, 3, 2, 0]
	npts_cdenom = 6

	# Faces: type, reference-to-face projection, normal

	faces = [

		('tri', lambda s, t: (t, s, -1), (0, 0, -1)),
		('tri', lambda s, t: ((t - 1)/2, s + (t + 1)/2, t), (-1, 0, 0.5)),
		('tri', lambda s, t: (s + (t + 1)/2, (t - 1)/2, t), (0, -1, 0.5)),
		('tri', lambda s, t: (s + (t + 1)/2, -s - (t + 1)/2, t), (1, 1, 0))
	]

	def __init__(self, nspts, cfg):
		super().__init__(nspts, cfg)
		# Move upper sol point to (0,0,1) (instead of (-1, -1, 1))
		rfac = 0.5*(self.upts[:,2] + 1.)
		self.upts[:,0] += rfac
		self.upts[:,1] += rfac
		# Flip upts/fpts
		self.upts[:,:2] *= -1
		self.fpts[:,:2] *= -1
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
		for i in range(self.nquadpts):
			coord = self.fpts[n, :]
			iidx = i % (self.order+1)
			negiidx = (self.order) - iidx
			jidx = i // (self.order+1)

			if negiidx == jidx:
				tag = 'm'
				idxl = np.where((tetl_bot == coord).all(axis=1))[0][0]
				idxr = np.where((tetr_bot == coord).all(axis=1))[0][0]
				assert idxl == idxr, f'{idxl}, {idxr}'
				self.fpts_idxs[n] = idxl
			elif negiidx > jidx:
				tag = 'l'
				self.fpts_idxs[n] = np.where((tetl_bot == coord).all(axis=1))[0][0]
			else:
				tag = 'r'
				self.fpts_idxs[n] = np.where((tetr_bot == coord).all(axis=1))[0][0]

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
			1: east
			2: north
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
				offset = 1*self.ntripts
				tet_coords = self.tetr.fpts[offset:offset+self.ntripts,:]
			elif face == 3:
				tag = 'r'
				offset = 2*self.ntripts
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
		newrow = np.zeros(self.nfpts)
		#Arbitrarily choose right side to include midpoints as well
		for i,v in enumerate(row):
			newidx = None
			if side == 'right':
				for j, vv in enumerate(self.fpts_idxs):
					if vv == i and self.fpts_map[j] in ['r', 'm']:
						newidx = j
			elif side == 'left':
				for j, vv in enumerate(self.fpts_idxs):
					if vv == i and self.fpts_map[j] in ['l', 'm']:
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
	def m3(self):
		M = np.zeros((self.nupts, self.nfpts))
		tol = 1e-8

		# Remove contribution of interior interface (tet diagonal)
		def remove_diag(m):
			return m[:-self.ntripts]


		for i, (x,y,z) in enumerate(self.upts):
			ml = remove_diag(self.tetl.gbasis_at([[x,y,z]])[0])
			mr = remove_diag(self.tetr.gbasis_at([[x,y,z]])[0])

			if np.abs(x + y) < tol: # Midpoints
				M[i, :] = 0.5*(self.expand_fpts_row(ml, 'left') + self.expand_fpts_row(mr, 'right'))
			elif (x + y) > tol: # Right
				M[i, :] = self.expand_fpts_row(mr, 'right')
			else: # Left 
				M[i, :] = self.expand_fpts_row(ml, 'left')
		return M

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
	pyr = PyrShape(False, cfg)
	m3 = pyr.m3
	tol = 1e-6

	u = np.zeros(pyr.nfpts)
	u[0] = 1.
	print(m3 @ u)


''' 
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