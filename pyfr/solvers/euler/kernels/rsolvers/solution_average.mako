# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures

% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(ul[{i}] + ur[{i}])'
                                 .format(i=i, j=j) for j in range(ndims))});
% endfor
</%pyfr:macro>
