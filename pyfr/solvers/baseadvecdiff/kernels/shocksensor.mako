# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='shocksensor' ndim='1'
              du='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              artvisc='out fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'>

// Calculate average solution defect
fpdtype_t int_du[${nvars}];
fpdtype_t max_du[${nvars}];
% for j in range(nvars):
    int_du[${j}] = 0.0;
    max_du[${j}] = 0.0;
% endfor
% for i,j in pyfr.ndrange(nupts, nvars):
    int_du[${j}] += ${weights[i]}*fabs(du[${i}][${j}]);
    max_du[${j}] = max(max_du[${j}], fabs(du[${i}][${j}]));
% endfor

// Calculate grid size
fpdtype_t h = 0.0;
% for i in range(nupts):
    h += ${weights[i]}*(1.0/rcpdjac[${i}]);
% endfor
h = pow(h, ${1.0/ndims})/${order + 1};

% for i,j in pyfr.ndrange(nupts, nvars):
    % if vis_method == 'pointwise':
        artvisc[${i}][${j}] = ${c_mu*vis_coeffs[j]}*pow(abs(du[${i}][${j}]), ${exp_fac})*h*h*${1.0/dt_rev};
    % elif vis_method == 'max':
        artvisc[${i}][${j}] = ${c_mu*vis_coeffs[j]}*pow(max_du[${j}], ${exp_fac})*h*h*${1.0/dt_rev};
    % elif vis_method == 'mean':
        artvisc[${i}][${j}] = ${c_mu*vis_coeffs[j]}*pow(int_du[${j}], ${exp_fac})*h*h*${1.0/dt_rev};
    % elif vis_method == 'constant':
        artvisc[${i}][${j}] = ${vis_coeffs[j]*c['mu_max']};
    % endif
    artvisc[${i}][${j}] = artvisc[${i}][${j}] < ${cutoff} ? 0.0 : artvisc[${i}][${j}];
% endfor

</%pyfr:kernel>
