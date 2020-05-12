# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.rsolvers.lambda'/>

<%pyfr:kernel name='shocksensor' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uf='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              artvisc='out fpdtype_t[${str(nupts)}][${str(nvars)}]'
              urcpdjac='in fpdtype_t[${str(nupts)}]'
              frcpdjac='in fpdtype_t[${str(nfpts)}]'
              usmats='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              rcpfsmats='in fpdtype_t[${str(nfpts)}][${str(ndims*ndims)}]'
              >

fpdtype_t ui[${nvars}], uj[${nvars}], uti[${nvars}], utj[${nvars}]; // Current point, adjacent point, transformed current point, transformed adjacent point
fpdtype_t c[${ndims}], cmag, cdelta[${ndims}], s; // C direction, c magnitude, surface c direction, max wave speed
fpdtype_t j11, j12, j13, j21, j22, j23; // Jacobian cofactor metrics (dxi/dx, etc.)
fpdtype_t fnorm[${ndims}], fmag; // Face normal, normal magnitude
fpdtype_t Jsurf; // Flux surface Jacobian
fpdtype_t djac;

% for i in range(nupts):
	% for k in range(nvars):
		artvisc[${i}][${k}] = 0.0;
	% endfor
	% for j in range(6):
		// If adjacent point is a flux point
		% if adjmat[i][j] >= nupts:  
			djac = 1.0/frcpdjac[${adjmat[i][j] - nupts}];
			% for k in range(nvars):
				ui[${k}] = u[${i}][${k}];
				uj[${k}] = uf[${adjmat[i][j] - nupts}][${k}];
			% endfor
			// If face is xi == += 1, use cofactors of eta,zeta. j1_ = eta, j2_ = zeta
			% if ffaces[adjmat[i][j] - nupts] == 2 or ffaces[adjmat[i][j] - nupts] == 4:
				j11 = djac*rcpfsmats[${adjmat[i][j] - nupts}][3]; j12 = djac*rcpfsmats[${adjmat[i][j] - nupts}][4]; j13 = djac*rcpfsmats[${adjmat[i][j] - nupts}][5]; 
				j21 = djac*rcpfsmats[${adjmat[i][j] - nupts}][6]; j22 = djac*rcpfsmats[${adjmat[i][j] - nupts}][7]; j23 = djac*rcpfsmats[${adjmat[i][j] - nupts}][8];
				% if ffaces[adjmat[i][j] - nupts] == 2: # xi = +1, fnorm = dx/dxi, dy/dxi, dz/dxi
					fnorm[0] = djac*rcpfsmats[${adjmat[i][j] - nupts}][0]; 
					fnorm[1] = djac*rcpfsmats[${adjmat[i][j] - nupts}][1]; 
					fnorm[2] = djac*rcpfsmats[${adjmat[i][j] - nupts}][2];
				% else: # xi = -1, fnorm = -dx/dxi, -dy/dxi, -dz/dxi
					fnorm[0] = -djac*rcpfsmats[${adjmat[i][j] - nupts}][0]; 
					fnorm[1] = -djac*rcpfsmats[${adjmat[i][j] - nupts}][1]; 
					fnorm[2] = -djac*rcpfsmats[${adjmat[i][j] - nupts}][2];
				% endif

			// If face is eta == += 1, use cofactors of xi,zeta. j1_ = xi, j2_ = zeta
			% elif ffaces[adjmat[i][j] - nupts] == 3 or ffaces[adjmat[i][j] - nupts] == 1:
				j11 = djac*rcpfsmats[${adjmat[i][j] - nupts}][0]; j12 = djac*rcpfsmats[${adjmat[i][j] - nupts}][1]; j13 = djac*rcpfsmats[${adjmat[i][j] - nupts}][2]; 
				j21 = djac*rcpfsmats[${adjmat[i][j] - nupts}][6]; j22 = djac*rcpfsmats[${adjmat[i][j] - nupts}][7]; j23 = djac*rcpfsmats[${adjmat[i][j] - nupts}][8];
				% if ffaces[adjmat[i][j] - nupts] == 3: # eta = +1, fnorm = dx/deta, dy/deta, dz/eta
					fnorm[0] = djac*rcpfsmats[${adjmat[i][j] - nupts}][3]; 
					fnorm[1] = djac*rcpfsmats[${adjmat[i][j] - nupts}][4]; 
					fnorm[2] = djac*rcpfsmats[${adjmat[i][j] - nupts}][5];
				% else: # eta = -1, fnorm = -dx/deta, -dy/deta, -dz/deta
					fnorm[0] = -djac*rcpfsmats[${adjmat[i][j] - nupts}][3]; 
					fnorm[1] = -djac*rcpfsmats[${adjmat[i][j] - nupts}][4]; 
					fnorm[2] = -djac*rcpfsmats[${adjmat[i][j] - nupts}][5];
				% endif

			// If face is zeta == += 1, use cofactors of xi, eta. j1_ = xi, j2_ = eta
			% elif ffaces[adjmat[i][j] - nupts] == 0 or ffaces[adjmat[i][j] - nupts] == 5:
				j11 = djac*rcpfsmats[${adjmat[i][j] - nupts}][0]; j12 = djac*rcpfsmats[${adjmat[i][j] - nupts}][1]; j13 = djac*rcpfsmats[${adjmat[i][j] - nupts}][2]; 
				j21 = djac*rcpfsmats[${adjmat[i][j] - nupts}][3]; j22 = djac*rcpfsmats[${adjmat[i][j] - nupts}][4]; j23 = djac*rcpfsmats[${adjmat[i][j] - nupts}][5];
				% if ffaces[adjmat[i][j] - nupts] == 2: # zeta = +1, fnorm = dx/dzeta, dy/dzeta, dz/dzeta
					fnorm[0] = djac*rcpfsmats[${adjmat[i][j] - nupts}][6]; 
					fnorm[1] = djac*rcpfsmats[${adjmat[i][j] - nupts}][7]; 
					fnorm[2] = djac*rcpfsmats[${adjmat[i][j] - nupts}][8];
				% else: # zeta = -1, fnorm = -dx/dzeta, -dy/dzeta, -dz/dzeta
					fnorm[0] = -djac*rcpfsmats[${adjmat[i][j] - nupts}][6]; 
					fnorm[1] = -djac*rcpfsmats[${adjmat[i][j] - nupts}][7]; 
					fnorm[2] = -djac*rcpfsmats[${adjmat[i][j] - nupts}][8];
				% endif
			% endif

			// Normalize face normal vectors
			fmag = 0;
			% for k in range(ndims):
				fmag += pow(fnorm[${k}], 2);
			% endfor
			fmag = pow(fmag, 0.5);

			% for k in range(ndims):
				fnorm[${k}] /= fmag;
			% endfor

			Jsurf = pow(j12*j23 - j13*j22, 2.0) + pow(j11*j23 - j13*j21, 2.0) + pow(j11*j22 - j12*j21, 2.0);
			Jsurf = pow(Jsurf, 0.5);

			//printf("%f\n", 55.388);
			//printf("%f\n", djac*rcpfsmats[${adjmat[i][j] - nupts}][0]);
			//printf("%f\n", djac*rcpfsmats[${adjmat[i][j] - nupts}][1]);
			//printf("%f\n", djac*rcpfsmats[${adjmat[i][j] - nupts}][2]);
			//printf("%f\n", djac*rcpfsmats[${adjmat[i][j] - nupts}][3]);
			//printf("%f\n", djac*rcpfsmats[${adjmat[i][j] - nupts}][4]);
			//printf("%f\n", djac*rcpfsmats[${adjmat[i][j] - nupts}][5]);
			//printf("%f\n", djac*rcpfsmats[${adjmat[i][j] - nupts}][6]);
			//printf("%f\n", djac*rcpfsmats[${adjmat[i][j] - nupts}][7]);
			//printf("%f\n", djac*rcpfsmats[${adjmat[i][j] - nupts}][8]);
			//printf("%f\n", 99.542);
			//printf("%f\n", frcpdjac[${adjmat[i][j] - nupts}]);



			% for k in range(ndims):
				cdelta[${k}] = ${prefac}*Jsurf*fnorm[${k}]/urcpdjac[${i}]; //*urcpdjac[${i}] ?
			% endfor


		// If adjacent point is a solution point
		% else:
			% for k in range(nvars):
				ui[${k}] = u[${i}][${k}];
				uj[${k}] = u[${adjmat[i][j]}][${k}];	
			% endfor
			% for k in range(ndims):
				cdelta[${k}] = 0.0;
			% endfor	
		% endif

		cmag = 0;
		% for k in range(ndims):
			c[${k}] = ${gradmat[i][j][0]}*usmats[${i}][${k}] + ${gradmat[i][j][1]}*usmats[${i}][${k+3}] + ${gradmat[i][j][2]}*usmats[${i}][${k+6}];
			c[${k}] /= urcpdjac[${i}]; // Because smats has a 1/det(J) term in it
			//c[${k}] -= cdelta[${k}];
			cmag += pow(c[${k}], 2);
		% endfor

		cmag = pow(cmag, 0.5);
		% for k in range(ndims):
			c[${k}] /= cmag; 
		% endfor

		${pyfr.expand('transform_to','c', 'ui', 'uti')}
    	${pyfr.expand('transform_to','c', 'uj', 'utj')}
    	${pyfr.expand('lambda_max','uti', 'utj', 's')}


		% for k in range(nvars):
			artvisc[${i}][${k}] += s*(uj[${k}] - ui[${k}])*cmag; 	
		% endfor
    %endfor

    %for k in range(nvars):
		//printf("%f\n", artvisc[${i}][${k}]);
	% endfor

%endfor
</%pyfr:kernel>
