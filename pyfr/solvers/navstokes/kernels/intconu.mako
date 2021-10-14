# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='intconu' ndim='1'
              ulin='in view fpdtype_t[${str(nvars)}]'
              urin='in view fpdtype_t[${str(nvars)}]'
              ulout='out view fpdtype_t[${str(nvars)}]'
              urout='out view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'>
% for i in range(nvars):
% if c['ldg-beta'] == -0.5:
    urout[${i}] = ulin[${i}];
% elif c['ldg-beta'] == 0.5:
    ulout[${i}] = urin[${i}];
% else:
    ulout[${i}] = urout[${i}] = urin[${i}]*${0.5 + c['ldg-beta']}
                              + ulin[${i}]*${0.5 - c['ldg-beta']};
% endif
% endfor

% if rsolver == 'upwind':
    fpdtype_t nv = ${pyfr.dot('nl[{i}]', 'ulin[{i}] + urin[{i}]', i=ndims)};
    fpdtype_t p = nv >= 0 ? ulin[${ndims}] : urin[${ndims}];
% else:
    fpdtype_t p;
    p = 0.5*(ulin[${ndims}] + urin[${ndims}]);
% endif

ulin[${ndims}] = p;
urin[${ndims}] = p;
ulout[${ndims}] = p;
urout[${ndims}] = p;

</%pyfr:kernel>
