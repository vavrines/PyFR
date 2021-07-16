# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='shocksensor' ndim='1'
              du='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              revvisc='out fpdtype_t[${str(nvars)}]'>

fpdtype_t int_du[${nvars}];
% for i,j in pyfr.ndrange(nupts, nvars):
    int_du[${j}] += ${weights[i]}*abs(du[${i}][${j}]);
% endfor

fpdtype_t mu = ${c['mu_max']};

% for i in range(nvars):
    revvisc[${i}] = mu;
% endfor


</%pyfr:kernel>
