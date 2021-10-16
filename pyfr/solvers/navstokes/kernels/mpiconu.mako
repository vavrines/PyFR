# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='mpiconu' ndim='1'
              ulin='inout view fpdtype_t[${str(nvars)}]'
              urin='inout mpi fpdtype_t[${str(nvars)}]'
              ulout='out view fpdtype_t[${str(nvars)}]'>
% for i in range(nvars):
% if c['ldg-beta'] == -0.5:
    ulout[${i}] = ulin[${i}];
% elif c['ldg-beta'] == 0.5:
    ulout[${i}] = urin[${i}];
% else:
    ulout[${i}] = urin[${i}]*${0.5 + c['ldg-beta']}
                + ulin[${i}]*${0.5 - c['ldg-beta']};
% endif
% endfor

fpdtype_t p = 0.5*(ulin[${ndims}] + urin[${ndims}]);
ulin[${ndims}] = p;
urin[${ndims}] = p;
ulout[${ndims}] = p;

</%pyfr:kernel>

