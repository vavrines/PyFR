# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='rdshocksensor' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              shockcell='out fpdtype_t'
              divf_fr='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              divf_rd='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              usmats='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'>


% if shocksensor == 'off':
	shockcell = 0;
% elif shocksensor == 'modal':
	<% se0 = crd['modal-sensor-coeff'] %>

	fpdtype_t totEn = 1e-15, pnEn = 1e-15, tmp;

	% for i, deg in enumerate(ubdegs):
		tmp = ${' + '.join('{jx}*u[{j}][{svar}]'.format(j=j, jx=jx, svar=svar)
							for j, jx in enumerate(invvdm[i]) if jx != 0)};
		totEn += tmp*tmp;
		% if deg >= order:
			pnEn += tmp*tmp;
		% endif
	% endfor

    fpdtype_t se  = pnEn/totEn;
    shockcell = (se < ${se0}) ? 0 : 1;
% elif shocksensor == 'maxmodal':
	<% se0 = crd['modal-sensor-coeff'] %>

	fpdtype_t totEn = 1e-15, pnEn = 1e-15, tmp;

	% for i, deg in enumerate(ubdegs):
		tmp = ${' + '.join('{jx}*u[{j}][{svar}]'.format(j=j, jx=jx, svar=svar)
							for j, jx in enumerate(invvdm[i]) if jx != 0)};
		totEn += tmp*tmp;
		% if deg == max(ubdegs):
			pnEn = tmp*tmp;
		% endif
	% endfor

    fpdtype_t se  = pnEn/totEn;
    shockcell = (se < ${se0}) ? 0 : 1;
% elif shocksensor == 'l1modal':
	<% se0 = crd['modal-sensor-coeff'] %>

	fpdtype_t totEn = 1e-15, pnEn = 1e-15, tmp;

	% for i, deg in enumerate(l1degs):
		tmp = ${' + '.join('{jx}*u[{j}][{svar}]'.format(j=j, jx=jx, svar=svar)
							for j, jx in enumerate(invvdm[i]) if jx != 0)};
		totEn += tmp*tmp;
		% if deg == max(l1degs):
			pnEn = max(tmp*tmp, pnEn);
		% endif
	% endfor

    fpdtype_t se  = pnEn/totEn;
    shockcell = (se < ${se0}) ? 0 : 1;

% elif shocksensor == 'conv-l2':
	<%include file='pyfr.solvers.euler.kernels.flux'/>
	fpdtype_t f1t[${nupts}][${ndims}][${nvars}], f2t[${nupts}][${ndims}][${nvars}], f3t[${nupts}][${ndims}][${nvars}];
	fpdtype_t f1[${nupts}][${ndims}][${nvars}], f2[${nupts}][${ndims}][${nvars}], f3[${nupts}][${ndims}][${nvars}];
	fpdtype_t u1[${nupts}][${nvars}], u2[${nupts}][${nvars}], u3[${nupts}][${nvars}];
	fpdtype_t f1c[${nupts}][${ndims}][${nvars}], f2c[${nupts}][${ndims}][${nvars}], f3c[${nupts}][${ndims}][${nvars}];
	fpdtype_t divf1[${nupts}][${nvars}], divf2[${nupts}][${nvars}], divf3[${nupts}][${nvars}];
	fpdtype_t dd12[${nvars}], dd23[${nvars}];

	// Calculate projected solution 
	% for i in range(nupts):
		% for var in range(nvars):
			u1[${i}][${var}] = ${' + '.join('{jx}*u[{j}][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(proj1[i]) if jx != 0)};
			u2[${i}][${var}] = ${' + '.join('{jx}*u[{j}][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(proj2[i]) if jx != 0)};
			u3[${i}][${var}] = ${' + '.join('{jx}*u[{j}][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(proj3[i]) if jx != 0)};
		% endfor
	% endfor

    fpdtype_t ftemp[${ndims}][${nvars}], utemp[${nvars}];
    fpdtype_t p, v[${ndims}];

	% for i in range(nupts):
		% for j in range(nvars):
			utemp[${j}] = u1[${i}][${j}];
		% endfor

    	${pyfr.expand('inviscid_flux', 'utemp', 'ftemp', 'p', 'v')};

		% for j in range(nvars):
			% for k in range(ndims):
				f1t[${i}][${k}][${j}] = ftemp[${k}][${j}];
			% endfor
		% endfor
	% endfor

	% for i in range(nupts):
		% for j in range(nvars):
			utemp[${j}] = u2[${i}][${j}];
		% endfor

    	${pyfr.expand('inviscid_flux', 'utemp', 'ftemp', 'p', 'v')};

		% for j in range(nvars):
			% for k in range(ndims):
				f2t[${i}][${k}][${j}] = ftemp[${k}][${j}];
			% endfor
		% endfor
	% endfor

	% for i in range(nupts):
		% for j in range(nvars):
			utemp[${j}] = u3[${i}][${j}];
		% endfor

    	${pyfr.expand('inviscid_flux', 'utemp', 'ftemp', 'p', 'v')};

		% for j in range(nvars):
			% for k in range(ndims):
				f3t[${i}][${k}][${j}] = ftemp[${k}][${j}];
			% endfor
		% endfor
	% endfor



	// Calculate projected flux values
	% for i in range(nupts):
		% for var,k in pyfr.ndrange(nvars, ndims):
			f1[${i}][${k}][${var}] = ${' + '.join('{jx}*f1t[{j}][{k}][{var}]'.format(j=j, jx=jx, k=k, var=var) for j, jx in enumerate(proj1[i]) if jx != 0)};
			f2[${i}][${k}][${var}] = ${' + '.join('{jx}*f2t[{j}][{k}][{var}]'.format(j=j, jx=jx, k=k, var=var) for j, jx in enumerate(proj2[i]) if jx != 0)};
			f3[${i}][${k}][${var}] = ${' + '.join('{jx}*f3t[{j}][{k}][{var}]'.format(j=j, jx=jx, k=k, var=var) for j, jx in enumerate(proj3[i]) if jx != 0)};
		% endfor
	% endfor

	// Transform to computational space
	% for i,dim,var in pyfr.ndrange(nupts, ndims,nvars):
    	f1c[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*f1[{0}][{2}][{var}]'.format(i, ndims*dim+k, k, var=var) for k in range(ndims))};
    	f2c[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*f2[{0}][{2}][{var}]'.format(i, ndims*dim+k, k, var=var) for k in range(ndims))};
    	f3c[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*f3[{0}][{2}][{var}]'.format(i, ndims*dim+k, k, var=var) for k in range(ndims))};
	% endfor

	// Calculate divergence (hardcoded for ndims=2)
	% for i,var in pyfr.ndrange(nupts, nvars):
		divf1[${i}][${var}] =  ${' + '.join('{jx}*f1c[{j}][0][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffxi[i]) if jx != 0)};
		divf2[${i}][${var}] =  ${' + '.join('{jx}*f2c[{j}][0][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffxi[i]) if jx != 0)};
		divf3[${i}][${var}] =  ${' + '.join('{jx}*f3c[{j}][0][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffxi[i]) if jx != 0)};
		divf1[${i}][${var}] += ${' + '.join('{jx}*f1c[{j}][1][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffeta[i]) if jx != 0)};
		divf2[${i}][${var}] += ${' + '.join('{jx}*f2c[{j}][1][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffeta[i]) if jx != 0)};
		divf3[${i}][${var}] += ${' + '.join('{jx}*f3c[{j}][1][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffeta[i]) if jx != 0)};
	% endfor

	% for i,var in pyfr.ndrange(nupts, nvars):
		divf1[${i}][${var}] *= rcpdjac[${i}]; 
		divf2[${i}][${var}] *= rcpdjac[${i}]; 
		divf3[${i}][${var}] *= rcpdjac[${i}]; 
	% endfor

	// L2 norm
	% for j in range(nvars):
		dd12[${j}] = 0;
		dd23[${j}] = 0;
	% endfor
	% for i,var in pyfr.ndrange(nupts, nvars):
		dd12[${var}] += ${quadwts[i]}*pow(divf1[${i}][${var}] - divf3[${i}][${var}], 2.0);
		dd23[${var}] += ${quadwts[i]}*pow(divf2[${i}][${var}] - divf3[${i}][${var}], 2.0);
	% endfor

	// Calculate characteristic length
	fpdtype_t h = 0.0;
	% for i in range(nupts):
		h += ${quadwts[i]}/rcpdjac[${i}];
	% endfor
	h = pow(h, 1.0/${ndims});

	shockcell = 0;
	// Look at density and energy
	% for i in [0,3]: 
		if (dd23[${i}] > (h*dd12[${i}] + ${crd['convergence-tolerance']})) { 
			shockcell = 1;
		}
	% endfor

% elif shocksensor == 'conv-maxdiv':
	<%include file='pyfr.solvers.euler.kernels.flux'/>
	fpdtype_t f1t[${nupts}][${ndims}][${nvars}], f2t[${nupts}][${ndims}][${nvars}], f3t[${nupts}][${ndims}][${nvars}];
	fpdtype_t f1[${nupts}][${ndims}][${nvars}], f2[${nupts}][${ndims}][${nvars}], f3[${nupts}][${ndims}][${nvars}];
	fpdtype_t u1[${nupts}][${nvars}], u2[${nupts}][${nvars}], u3[${nupts}][${nvars}];
	fpdtype_t f1c[${nupts}][${ndims}][${nvars}], f2c[${nupts}][${ndims}][${nvars}], f3c[${nupts}][${ndims}][${nvars}];
	fpdtype_t divf1[${nupts}][${nvars}], divf2[${nupts}][${nvars}], divf3[${nupts}][${nvars}];
	fpdtype_t dd12[${nvars}], dd23[${nvars}];

	// Calculate projected solution 
	% for i in range(nupts):
		% for var in range(nvars):
			u1[${i}][${var}] = ${' + '.join('{jx}*u[{j}][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(proj1[i]) if jx != 0)};
			u2[${i}][${var}] = ${' + '.join('{jx}*u[{j}][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(proj2[i]) if jx != 0)};
			u3[${i}][${var}] = ${' + '.join('{jx}*u[{j}][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(proj3[i]) if jx != 0)};
		% endfor
	% endfor

    fpdtype_t ftemp[${ndims}][${nvars}], utemp[${nvars}];
    fpdtype_t p, v[${ndims}];

	% for i in range(nupts):
		% for j in range(nvars):
			utemp[${j}] = u1[${i}][${j}];
		% endfor

    	${pyfr.expand('inviscid_flux', 'utemp', 'ftemp', 'p', 'v')};

		% for j in range(nvars):
			% for k in range(ndims):
				f1t[${i}][${k}][${j}] = ftemp[${k}][${j}];
			% endfor
		% endfor
	% endfor

	% for i in range(nupts):
		% for j in range(nvars):
			utemp[${j}] = u2[${i}][${j}];
		% endfor

    	${pyfr.expand('inviscid_flux', 'utemp', 'ftemp', 'p', 'v')};

		% for j in range(nvars):
			% for k in range(ndims):
				f2t[${i}][${k}][${j}] = ftemp[${k}][${j}];
			% endfor
		% endfor
	% endfor

	% for i in range(nupts):
		% for j in range(nvars):
			utemp[${j}] = u3[${i}][${j}];
		% endfor

    	${pyfr.expand('inviscid_flux', 'utemp', 'ftemp', 'p', 'v')};

		% for j in range(nvars):
			% for k in range(ndims):
				f3t[${i}][${k}][${j}] = ftemp[${k}][${j}];
			% endfor
		% endfor
	% endfor



	// Calculate projected flux values
	% for i in range(nupts):
		% for var,k in pyfr.ndrange(nvars, ndims):
			f1[${i}][${k}][${var}] = ${' + '.join('{jx}*f1t[{j}][{k}][{var}]'.format(j=j, jx=jx, k=k, var=var) for j, jx in enumerate(proj1[i]) if jx != 0)};
			f2[${i}][${k}][${var}] = ${' + '.join('{jx}*f2t[{j}][{k}][{var}]'.format(j=j, jx=jx, k=k, var=var) for j, jx in enumerate(proj2[i]) if jx != 0)};
			f3[${i}][${k}][${var}] = ${' + '.join('{jx}*f3t[{j}][{k}][{var}]'.format(j=j, jx=jx, k=k, var=var) for j, jx in enumerate(proj3[i]) if jx != 0)};
		% endfor
	% endfor

	// Transform to computational space
	% for i,dim,var in pyfr.ndrange(nupts, ndims,nvars):
    	f1c[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*f1[{0}][{2}][{var}]'.format(i, ndims*dim+k, k, var=var) for k in range(ndims))};
    	f2c[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*f2[{0}][{2}][{var}]'.format(i, ndims*dim+k, k, var=var) for k in range(ndims))};
    	f3c[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*f3[{0}][{2}][{var}]'.format(i, ndims*dim+k, k, var=var) for k in range(ndims))};
	% endfor

	// Set max div to 0
	fpdtype_t maxd1[${nvars}], maxd2[${nvars}], maxd3[${nvars}];
	% for var in range(nvars):
		maxd1[${var}] = 0.0;
		maxd2[${var}] = 0.0;
		maxd3[${var}] = 0.0;
	% endfor

	// Calculate divergence (hardcoded for ndims=2)
	fpdtype_t tmp1, tmp2;
	% for i,var in pyfr.ndrange(nupts, nvars):
		tmp1 = rcpdjac[${i}]*(${' + '.join('{jx}*f1c[{j}][0][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffxi[i]) if jx != 0)});
		tmp2 = rcpdjac[${i}]*(${' + '.join('{jx}*f1c[{j}][1][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffeta[i]) if jx != 0)});
		maxd1[${var}] = max(maxd1[${var}], pow(tmp1*tmp1 + tmp2*tmp2, 2.0));
		tmp1 = rcpdjac[${i}]*(${' + '.join('{jx}*f2c[{j}][0][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffxi[i]) if jx != 0)});
		tmp2 = rcpdjac[${i}]*(${' + '.join('{jx}*f2c[{j}][1][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffeta[i]) if jx != 0)});		
		maxd2[${var}] = max(maxd2[${var}], pow(tmp1*tmp1 + tmp2*tmp2, 2.0));
		tmp1 = rcpdjac[${i}]*(${' + '.join('{jx}*f3c[{j}][0][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffxi[i]) if jx != 0)});
		tmp2 = rcpdjac[${i}]*(${' + '.join('{jx}*f3c[{j}][1][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffeta[i]) if jx != 0)});
		maxd3[${var}] = max(maxd3[${var}], pow(tmp1*tmp1 + tmp2*tmp2, 2.0));
	% endfor

	shockcell = 0;
	// Look at density and energy
	% for i in [0,1,2]: 
		if (maxd2[${i}] > (maxd1[${i}] + ${crd['convergence-tolerance']})  &&
			maxd3[${i}] > (maxd2[${i}] + ${crd['convergence-tolerance']})) { 
				shockcell = 1;
		}
	% endfor

	
% elif shocksensor == 'conv-mag':
	<%include file='pyfr.solvers.euler.kernels.flux'/>
	fpdtype_t f1t[${nupts}][${ndims}][${nvars}], f2t[${nupts}][${ndims}][${nvars}], f3t[${nupts}][${ndims}][${nvars}];
	fpdtype_t f1[${nupts}][${ndims}][${nvars}], f2[${nupts}][${ndims}][${nvars}], f3[${nupts}][${ndims}][${nvars}];
	fpdtype_t u1[${nupts}][${nvars}], u2[${nupts}][${nvars}], u3[${nupts}][${nvars}];
	fpdtype_t f1c[${nupts}][${ndims}][${nvars}], f2c[${nupts}][${ndims}][${nvars}], f3c[${nupts}][${ndims}][${nvars}];
	fpdtype_t divf1[${nupts}][${nvars}], divf2[${nupts}][${nvars}], divf3[${nupts}][${nvars}];
	fpdtype_t dd12[${nvars}], dd23[${nvars}];

	// Calculate projected solution 
	% for i in range(nupts):
		% for var in range(nvars):
			u1[${i}][${var}] = ${' + '.join('{jx}*u[{j}][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(proj1[i]) if jx != 0)};
			u2[${i}][${var}] = ${' + '.join('{jx}*u[{j}][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(proj2[i]) if jx != 0)};
			u3[${i}][${var}] = ${' + '.join('{jx}*u[{j}][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(proj3[i]) if jx != 0)};
		% endfor
	% endfor

    fpdtype_t ftemp[${ndims}][${nvars}], utemp[${nvars}];
    fpdtype_t p, v[${ndims}];

	% for i in range(nupts):
		% for j in range(nvars):
			utemp[${j}] = u1[${i}][${j}];
		% endfor

    	${pyfr.expand('inviscid_flux', 'utemp', 'ftemp', 'p', 'v')};

		% for j in range(nvars):
			% for k in range(ndims):
				f1t[${i}][${k}][${j}] = ftemp[${k}][${j}];
			% endfor
		% endfor
	% endfor

	% for i in range(nupts):
		% for j in range(nvars):
			utemp[${j}] = u2[${i}][${j}];
		% endfor

    	${pyfr.expand('inviscid_flux', 'utemp', 'ftemp', 'p', 'v')};

		% for j in range(nvars):
			% for k in range(ndims):
				f2t[${i}][${k}][${j}] = ftemp[${k}][${j}];
			% endfor
		% endfor
	% endfor

	% for i in range(nupts):
		% for j in range(nvars):
			utemp[${j}] = u3[${i}][${j}];
		% endfor

    	${pyfr.expand('inviscid_flux', 'utemp', 'ftemp', 'p', 'v')};

		% for j in range(nvars):
			% for k in range(ndims):
				f3t[${i}][${k}][${j}] = ftemp[${k}][${j}];
			% endfor
		% endfor
	% endfor



	// Calculate projected flux values
	% for i in range(nupts):
		% for var,k in pyfr.ndrange(nvars, ndims):
			f1[${i}][${k}][${var}] = ${' + '.join('{jx}*f1t[{j}][{k}][{var}]'.format(j=j, jx=jx, k=k, var=var) for j, jx in enumerate(proj1[i]) if jx != 0)};
			f2[${i}][${k}][${var}] = ${' + '.join('{jx}*f2t[{j}][{k}][{var}]'.format(j=j, jx=jx, k=k, var=var) for j, jx in enumerate(proj2[i]) if jx != 0)};
			f3[${i}][${k}][${var}] = ${' + '.join('{jx}*f3t[{j}][{k}][{var}]'.format(j=j, jx=jx, k=k, var=var) for j, jx in enumerate(proj3[i]) if jx != 0)};
		% endfor
	% endfor

	// Transform to computational space
	% for i,dim,var in pyfr.ndrange(nupts, ndims,nvars):
    	f1c[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*f1[{0}][{2}][{var}]'.format(i, ndims*dim+k, k, var=var) for k in range(ndims))};
    	f2c[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*f2[{0}][{2}][{var}]'.format(i, ndims*dim+k, k, var=var) for k in range(ndims))};
    	f3c[${i}][${dim}][${var}] = ${' + '.join('usmats[{0}][{1}]*f3[{0}][{2}][{var}]'.format(i, ndims*dim+k, k, var=var) for k in range(ndims))};
	% endfor

	// Calculate divergence (hardcoded for ndims=2)
	% for i,var in pyfr.ndrange(nupts, nvars):
		divf1[${i}][${var}] =  pow((${' + '.join('{jx}*f1c[{j}][0][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffxi[i]) if jx != 0)}), 2.0)
								+ pow((${' + '.join('{jx}*f1c[{j}][1][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffeta[i]) if jx != 0)}), 2.0);
		divf2[${i}][${var}] =  pow((${' + '.join('{jx}*f2c[{j}][0][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffxi[i]) if jx != 0)}), 2.0)
								+ pow((${' + '.join('{jx}*f2c[{j}][1][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffeta[i]) if jx != 0)}), 2.0);
		divf3[${i}][${var}] =  pow((${' + '.join('{jx}*f3c[{j}][0][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffxi[i]) if jx != 0)}), 2.0)
								+ pow((${' + '.join('{jx}*f3c[{j}][1][{var}]'.format(j=j, jx=jx, var=var) for j, jx in enumerate(diffeta[i]) if jx != 0)}), 2.0);
		divf1[${i}][${var}] = pow(divf1[${i}][${var}], 0.5);
		divf2[${i}][${var}] = pow(divf2[${i}][${var}], 0.5);
		divf3[${i}][${var}] = pow(divf3[${i}][${var}], 0.5);
	% endfor

	// L2 norm
	% for j in range(nvars):
		dd12[${j}] = 0;
		dd23[${j}] = 0;
	% endfor
	% for i,var in pyfr.ndrange(nupts, nvars):
		dd12[${var}] += ${quadwts[i]}*pow(divf1[${i}][${var}] - divf2[${i}][${var}], 2.0);
		dd23[${var}] += ${quadwts[i]}*pow(divf2[${i}][${var}] - divf3[${i}][${var}], 2.0);
	% endfor

	// Calculate characteristic length
	fpdtype_t h = 0.0;
	% for i in range(nupts):
		h += ${quadwts[i]}/rcpdjac[${i}];
	% endfor
	h = pow(h, 1.0/${ndims});

	shockcell = 0;
	// Look at density and energy
	% for i in [0,3]: 
		if (dd23[${i}] > (h*dd12[${i}] + ${crd['convergence-tolerance']})) { 
			shockcell = 1;
		}
	% endfor

% elif shocksensor == 'gradient':
	fpdtype_t gradrho[${nupts}];
	% for i in range(nupts):
		gradrho[${i}] =  pow(${' + '.join('{jx}*u[{j}][0]'.format(j=j, jx=jx) for j, jx in enumerate(diffxi[i]) if jx != 0)},  2.0);
		gradrho[${i}] += pow(${' + '.join('{jx}*u[{j}][0]'.format(j=j, jx=jx) for j, jx in enumerate(diffeta[i]) if jx != 0)}, 2.0);
		gradrho[${i}] = pow(gradrho[${i}], 0.5);
	% endfor

	shockcell = 0;
	% for i in range(nupts): 
		if (gradrho[${i}] > ${crd['max-grad']}) { 
			shockcell = 1;
		}
	% endfor

% elif shocksensor == 'averagegradient':
	fpdtype_t gradrho = 0, tmp;
	% for i in range(nupts):
		tmp =  pow(${' + '.join('{jx}*u[{j}][0]'.format(j=j, jx=jx) for j, jx in enumerate(diffxi[i]) if jx != 0)},  2.0);
		tmp += pow(${' + '.join('{jx}*u[{j}][0]'.format(j=j, jx=jx) for j, jx in enumerate(diffeta[i]) if jx != 0)}, 2.0);
		gradrho += pow(tmp, 0.5)/${nupts};
	% endfor

	shockcell = 0;
	if (gradrho > ${crd['max-grad']} || gradrho < ${crd['min-grad']}) { 
		shockcell = 1;
	}

% else:
	shockcell = 1;
% endif

if (shockcell == 1) {
	% for i in range(nupts):
		% for j in range(nvars):
			divf_fr[${i}][${j}] = divf_rd[${i}][${j}];
		% endfor
	% endfor
}

</%pyfr:kernel>