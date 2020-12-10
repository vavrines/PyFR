# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='riemanndifference' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              plocu='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              usmats='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              uf='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              fsmats='in fpdtype_t[${str(nfpts)}][${str(ndims*ndims)}]'
              divf='out fpdtype_t[${str(nupts)}][${str(nvars)}]'
              res='in fpdtype_t[${str(nupts)}]'
              >

<%pyfr:macro name='get_normal' params='xl,xr,n'>
	fpdtype_t nmag = 0;
	% for dim in range(ndims):
		n[${dim}] = xl[${dim}] - xr[${dim}];
	% endfor		

	// Normalize norm
	% for dim in range(ndims):
		nmag += pow(n[${dim}], 2);
	% endfor
	nmag = pow(nmag, 0.5);
	% for dim in range(ndims):
		n[${dim}] /= nmag; 
	% endfor
</%pyfr:macro>

<%pyfr:macro name='get_tangent2d' params='xl,xr,t'>
	// Counter-clockwise rotation
	fpdtype_t tmp;
	${pyfr.expand('get_normal','xl', 'xr', 't')}

	tmp = t[1];
	t[1] = -t[0];
	t[0] = tmp;
</%pyfr:macro>

<%pyfr:macro name='centered_rsolve' params='ul, ur, t, f'>
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}], p, v[${ndims}];
    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'p', 'v')}
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'p', 'v')}

    % for var in range(nvars):
        f[${var}] = 0.5*(${' + '.join('t[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=var, j=j) for j in range(ndims))});
    % endfor
    
</%pyfr:macro>

% for i,j in pyfr.ndrange(nupts, nvars):
	divf[${i}][${j}] = 0.0;
% endfor


fpdtype_t ftemp[${ndims}][${nvars}], ftemp2[${ndims}][${nvars}], fntemp[${nvars}], fntemp2[${nvars}], xl[${ndims}], xr[${ndims}];
fpdtype_t line_sol[${order+1}][${nvars}], line_flux[${order+2}][${ndims}][${nvars}], line_tflux[${order+2}][${ndims}][${nvars}];
fpdtype_t ul[${nvars}], ur[${nvars}], usd[${nvars}], n[${ndims}], t[${ndims}], tmp, p, v[${ndims}];

// Perform along xi direction
% for j in range(order+1):
	// Gather solution along constant eta line
	% for i, var in pyfr.ndrange(order+1, nvars):
		line_sol[${i}][${var}] = u[${i+j*(order+1)}][${var}];
	% endfor

	% for i in range(order):
		// Get normal direction between solution points
		% for dim in range(ndims):
			xl[${dim}] = plocu[${i+1+j*(order+1)}][${dim}];
			xr[${dim}] = plocu[${i+  j*(order+1)}][${dim}];
		% endfor		
		${pyfr.expand('get_normal','xl', 'xr', 'n')}
		${pyfr.expand('get_tangent2d','xl', 'xr', 't')}

		// Set left and right solution states
		% for var in range(nvars):
			ul[${var}] = line_sol[${i}][${var}];
			ur[${var}] = line_sol[${i+1}][${var}];
		% endfor

        // Interpolate to RD points
        % for var in range(nvars):
            usd[${var}] = ${' + '.join('{mx}*line_sol[{m}][{var}]'.format(m=m, mx=mx, var=var)
                       for m, mx in enumerate(interpmat[i+1]) if mx != 0)}; // 
        % endfor

        // Check if density or pressure is negative 
        p = ${c['gamma'] - 1}*(usd[${nvars-2}] - 0.5*(pow(usd[1], 2.0) + pow(usd[2], 2.0))/usd[0]);

        // If res > e_max or rho < tol or pressure < tol
        if (abs(res[${i+j*(order+1)}]) > ${e_max} 
            || abs(res[${i+j*(order+1)+1}]) > ${e_max}
            || usd[0] < ${tol}
            || p < ${tol}) {

            ${pyfr.expand('rsolve','ul','ur','n','fntemp')};
            ${pyfr.expand('centered_rsolve','ul','ur','t','fntemp2')};
            % for var in range(nvars):
                line_flux[${i+1}][0][${var}] = n[0]*fntemp[${var}] + n[1]*fntemp2[${var}];
                line_flux[${i+1}][1][${var}] = n[1]*fntemp[${var}] - n[0]*fntemp2[${var}];
            % endfor
            }
        else {
            // And calculate flux
            ${pyfr.expand('inviscid_flux', 'usd', 'ftemp', 'p', 'v')};
            % for dim, var in pyfr.ndrange(ndims, nvars):
                line_flux[${i+1}][${dim}][${var}] = ftemp[${dim}][${var}];
            % endfor
        }
    % endfor

    // Transform flux to computational space
	// Definitely ndims*dim+k
	% for i in range(1,order+1):
		% for dim, var in pyfr.ndrange(ndims, nvars):
			line_tflux[${i}][${dim}][${var}] = ${' + '.join('(0.5*usmats[{0}][{2}] + 0.5*usmats[{1}][{2}])*line_flux[{5}][{3}][{4}]'
                                                 .format(i-1+j*(order+1),i+j*(order+1), ndims*dim+k, k, var, i) for k in range(ndims))};
		% endfor
	% endfor

    // Set interface fluxes to zero (include them in correction term)
    % for dim, var in pyfr.ndrange(ndims, nvars):
        line_tflux[0][${dim}][${var}] = 0;
        line_tflux[${order+1}][${dim}][${var}] = 0;
    % endfor

	// Calculate df/dxi at solution points
	% for var, i in pyfr.ndrange(nvars, order+1):
		tmp =  ${' + '.join('{mx}*line_tflux[{m}][0][{var}]'.format(m=m, mx=mx, var=var)
                   for m, mx in enumerate(diffmatRD[i]) if mx != 0)};
		divf[${i+ j*(order+1)}][${var}] += tmp;
	% endfor 
% endfor

// Perform along eta direction
% for i in range(order+1):
	// Gather solution along constant xi line
	% for j, var in pyfr.ndrange(order+1, nvars):
		line_sol[${j}][${var}] = u[${i+j*(order+1)}][${var}];
	% endfor

	% for j in range(order):
		// Get normal direction between solution points
		% for dim in range(ndims):
			xl[${dim}] = plocu[${i+(j+1)*(order+1)}][${dim}];
			xr[${dim}] = plocu[${i+j*(order+1)}][${dim}];
		% endfor		
		${pyfr.expand('get_normal','xl', 'xr', 'n')}
		${pyfr.expand('get_tangent2d','xl', 'xr', 't')}

		// Set left and right solution states
		% for var in range(nvars):
			ul[${var}] = line_sol[${j}][${var}];
			ur[${var}] = line_sol[${j+1}][${var}];
		% endfor

        // Interpolate to RD points
        % for var in range(nvars):
            usd[${var}] = ${' + '.join('{mx}*line_sol[{m}][{var}]'.format(m=m, mx=mx, var=var)
                       for m, mx in enumerate(interpmat[j+1]) if mx != 0)};
        % endfor

        // Check if density or pressure is negative 
        p = ${c['gamma'] - 1}*(usd[${nvars-2}] - 0.5*(pow(usd[1], 2.0) + pow(usd[2], 2.0))/usd[0]);

        // If res > e_max or rho < tol or pressure < tol
        if (abs(res[${i+(j+1)*(order+1)}]) > ${e_max} 
            || abs(res[${i+j*(order+1)}]) > ${e_max}
            || usd[0] < ${tol}
            || p < ${tol}) {

    		${pyfr.expand('rsolve','ul','ur','n','fntemp')};
            ${pyfr.expand('centered_rsolve','ul','ur','t','fntemp2')};
    		% for var in range(nvars):
    			line_flux[${j+1}][0][${var}] = n[0]*fntemp[${var}] + n[1]*fntemp2[${var}];
    			line_flux[${j+1}][1][${var}] = n[1]*fntemp[${var}] - n[0]*fntemp2[${var}];
    		% endfor
            }
        else {
            // And calculate flux
            ${pyfr.expand('inviscid_flux', 'usd', 'ftemp', 'p', 'v')};
            % for dim, var in pyfr.ndrange(ndims, nvars):
                line_flux[${j+1}][${dim}][${var}] = ftemp[${dim}][${var}];
            % endfor
        }
	% endfor

	// Transform flux to computational space
	% for j in range(1,order+1):
		% for dim, var in pyfr.ndrange(ndims, nvars):
			line_tflux[${j}][${dim}][${var}] = ${' + '.join('(0.5*usmats[{0}][{2}] + 0.5*usmats[{1}][{2}])*line_flux[{5}][{3}][{4}]'
                                                 .format(i+(j-1)*(order+1),i+j*(order+1), ndims*dim+k, k, var, j) for k in range(ndims))};
		% endfor
	% endfor

    // Set interface fluxes to zero (include them in correction term)
    % for dim, var in pyfr.ndrange(ndims, nvars):
        line_tflux[0][${dim}][${var}] = 0;
        line_tflux[${order+1}][${dim}][${var}] = 0;
    % endfor

	// Calculate dg/deta at solution points
	% for var, j in pyfr.ndrange(nvars, order+1):
		tmp =  ${' + '.join('{mx}*line_tflux[{m}][1][{var}]'.format(m=m, mx=mx, var=var)
                   for m, mx in enumerate(diffmatRD[j]) if mx != 0)};
		divf[${i+ j*(order+1)}][${var}] += tmp; 
	% endfor

% endfor





</%pyfr:kernel>
