# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<%pyfr:kernel name='shocksensor' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              artvisc='out fpdtype_t[${str(nvars)}]'>



fpdtype_t mu = ${c['max-artvisc']};

% for i in range(nvars):
    artvisc[${i}] = mu;
% endfor

</%pyfr:kernel>
