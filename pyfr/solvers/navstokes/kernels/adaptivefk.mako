# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='adaptivefk' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              fk='out fpdtype_t'>

// Evaluate modal coefficients
fpdtype_t k_res = ${c['min_ku']}, tmp;

% for i, deg in enumerate(ubdegs):
	% for k in range(1,ndims+1):
    	tmp = ${' + '.join('{jx}*u[{j}][{k}]'.format(j=j, jx=jx, k=k)
                       for j, jx in enumerate(invvdm[i]) if jx != 0)};
		% if deg > 0:
		    k_res += tmp*tmp;
		% endif
	% endfor
% endfor

fpdtype_t djac = 0;
fpdtype_t rho = 0;
fpdtype_t wu = 0;
fpdtype_t ku = 0;

// Average rho, djac, and wu over element
% for i in range(nupts):
	djac += (1.0/rcpdjac[${i}])/${nupts}; 
	rho  += u[${i}][0]/${nupts};
	wu   += exp(u[${i}][${nvars-1}])/${nupts};
	ku   += u[${i}][${nvars-2}]/${nupts};
% endfor


// Calculate minimum turbulent length scale = (nu^3/epsilon)^(1/4) = (nu^3/(betastar*k*w))^(1/4)
fpdtype_t dx_max = pow(${c['elemvol']}*djac, 0.33); 
dx_max = dx_max/(${order + 1}); // Divide by order+1 to account for increased resolving power at higher-order
fpdtype_t l_min   = pow(((pow(${c['mu']}/rho, 3))/(${c['betastar']}*k_res*wu)), 0.25); // Dzanic method
fpdtype_t l_kw    = pow(ku,0.5)/(${c['betastar']}*wu); // Girimaji method
fpdtype_t l_kwhyb = pow(k_res,0.5)/(${c['betastar']}*wu); // Hybrid method

fpdtype_t l_scale = 0;
% if adpans == 'dzanic':
	l_scale = l_min;
% elif adpans == 'girimaji':
	l_scale = l_kw;
% elif adpans == 'hybrid':
	l_scale = l_kwhyb;
% endif

// Calculate variable fk and assume fw = 1/fk
fpdtype_t fk_tar = ${c['C_PANS']}*pow(dx_max/l_min, 2.0/3.0);
fk = fk + ${c['fktdc']}*(fk_tar - fk);


</%pyfr:kernel>



