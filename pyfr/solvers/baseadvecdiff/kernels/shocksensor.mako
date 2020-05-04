# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.rsolvers.lambda'/>

<% se0 = math.log10(cav['s0']) %>

<%pyfr:kernel name='shocksensor' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uf='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              plocu='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              plocf='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              artvisc='out fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='out fpdtype_t[${str(nupts)}]'
              >
    // Smoothness indicator
    fpdtype_t totEn = 0.0, pnEn = 1e-15, tmp;

% for i, deg in enumerate(ubdegs):
    tmp = ${' + '.join('{jx}*u[{j}][{svar}]'.format(j=j, jx=jx, svar=svar)
                       for j, jx in enumerate(invvdm[i]) if jx != 0)};

    totEn += tmp*tmp;
% if deg >= order:
    pnEn += tmp*tmp;
% endif
% endfor

fpdtype_t se  = ${1/math.log(10)}*log(pnEn/totEn);

// Compute cell-wise artificial viscosity
//fpdtype_t mu = (se < ${se0 - cav['kappa']})
//             ? 0.0
//             : ${0.5*cav['max-artvisc']}*(1.0 + sin(${0.5*math.pi/cav['kappa']}*(se - ${se0})));
//mu = (se < ${se0 + cav['kappa']}) ? mu : ${cav['max-artvisc']};


fpdtype_t ui[${nvars}], uj[${nvars}], uti[${nvars}], utj[${nvars}]; // Current point, adjacent point, transformed current point, transformed adjacent point
fpdtype_t n[${ndims}], nmag, s; // Normal direction, magnitude, max wave speed

% for i in range(nupts):
	% for k in range(nvars):
		artvisc[${i}][${k}] = 0;
	% endfor
	% for j in range(6):
		// If adjacent point is a flux point
		% if adjmat[i][j] > nupts:  
			% for k in range(nvars):
				ui[${k}] = u[${i}][${k}];
				uj[${k}] = uf[${adjmat[i][j] - nupts}][${k}];
			% endfor
			% for k in range(ndims):
				n[${k}] = plocf[${adjmat[i][j] - nupts}][${k}] - plocu[${i}][${k}];
			% endfor
			//printf("%f\n", 888.8);
			//printf("%f\n", plocf[${adjmat[i][j] - nupts}][0]);
			//printf("%f\n", plocf[${adjmat[i][j] - nupts}][1]);
			//printf("%f\n", plocf[${adjmat[i][j] - nupts}][2]);


		// If adjacent point is a solution point
		% else:
			% for k in range(nvars):
				ui[${k}] = u[${i}][${k}];
				uj[${k}] = u[${adjmat[i][j]}][${k}];	
			% endfor		
			% for k in range(ndims):
				n[${k}] = plocu[${adjmat[i][j]}][${k}] - plocu[${i}][${k}];
			% endfor
			//printf("%f\n", 44.4);
			//printf("%f\n", plocu[${adjmat[i][j]}][0]);
			//printf("%f\n", plocu[${adjmat[i][j]}][1]);
			//printf("%f\n", plocu[${adjmat[i][j]}][2]);
		% endif
		//printf("%f\n", plocu[${i}][0]);
		//printf("%f\n", plocu[${i}][1]);
		//printf("%f\n", plocu[${i}][2]);
		//printf("%f\n", n[0]);
		//printf("%f\n", n[1]);
		//printf("%f\n", n[2]);
		//printf("%f\n", 69.9);

		// Normalize norm
		nmag = 0;
		% for k in range(ndims):
			nmag += pow(n[${k}], 2);
		% endfor
		nmag = pow(nmag, 0.5);
		% for k in range(ndims):
			n[${k}] /= nmag; 
		% endfor

		${pyfr.expand('transform_to','n', 'ui', 'uti')}
    	${pyfr.expand('transform_to','n', 'uj', 'utj')}
    	${pyfr.expand('lambda_max','uti', 'utj', 's')}

		% for k in range(nvars):
			artvisc[${i}][${k}] += ${weights[i]}*s*(uj[${k}] - ui[${k}])*nmag/rcpdjac[${i}];
		% endfor

    %endfor
    %for k in range(nvars):
		//printf("%f\n", artvisc[${i}][${k}]);
	% endfor

%endfor
</%pyfr:kernel>
