from pyfr.shapes import *
from pyfr.polys import *
from pyfr.inifile import Inifile



path = '/mnt/c/Users/Tarik/Documents/GitHub/PyFR-PANS/pyfr/test.ini'
cfg = Inifile.load(path)


class TetLShape(BaseShape):

	'''
	This tet defined by points
	
	(-1, -1)
	( 1, -1)
	(-1,  1)
	( 0,  0)

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
			pnew = (1 - rfac)*p + rfac*(0.0)
			qnew = (1 - rfac)*q + rfac*(0.0)
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
		self.upts[:,:2] *= -1
		self.fpts[:,:2] *= -1
		self.ubasis = get_polybasis(self.name, self.order + 1, self.upts)

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
			pnew = (1 - rfac)*p + rfac*(0.0)
			qnew = (1 - rfac)*q + rfac*(0.0) 
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
		print(self.tetr.fpts)
		print()
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



		print(self.fpts_map)
		print(self.fpts_idxs)

	def expand_upts_row(row, side):
		pass

	def expand_fpts_row(row, side):
		pass
	
	@cached_property
	def m0(self):
		M = np.zeros((self.nfpts, self.nupts))
		m0l = self.tetl.m0
		for i, (x,y,z) in enumerate(self.fpts):
			if (x + y) > 1e-8: # Strict-right
				row = self.tetr.ubasis.nodal_basis_at([[x,y,z]])
			else: # Left or diag
				row = self.tetl.ubasis.nodal_basis_at([[x,y,z]])
			# print(row)
			M[i, :] = row
		return M


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

# print(pyr.std_ele(1))
# print(tetl.std_ele(1))
# print(tetr.std_ele(1))

# print()
# print(tetl.fpts)
# print(tetr.fpts)
# print(pyr.fpts)
# print()
# print(pyr.m0)
