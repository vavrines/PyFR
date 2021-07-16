# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='shocksensor' ndim='1'
              du='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              revvisc='out fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'>

// Calculate average solution defect
fpdtype_t int_du[${nvars}];
% for j in range(nvars):
    int_du[${j}] = 0.0;
% endfor
% for i,j in pyfr.ndrange(nupts, nvars):
    int_du[${j}] += ${weights[i]}*fabs(du[${i}][${j}]);
% endfor

// Calculate grid size
fpdtype_t h = 0.0;
% for i in range(nupts):
    h += ${weights[i]}*(1.0/rcpdjac[${i}]);
% endfor
h = pow(h, ${1.0/ndims})/${order + 1};


% for i,j in pyfr.ndrange(nupts, nvars):
    revvisc[${i}][${j}] = fmin(${c_mu*vis_coeffs[j]}*int_du[${j}]*h*h*${1.0/dt_rev}, ${c['mu_max']});
% endfor


</%pyfr:kernel>
