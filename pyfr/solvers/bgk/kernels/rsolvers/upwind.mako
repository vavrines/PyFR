# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='rsolve' params='fl, fr, n, nF'>
    fpdtype_t nFl, nFr, nv;

% for j in range(nvars):
    nFl = nFr = nv = 0.0;
    % for i in range(ndims):
    nFl += n[${i}]*(${-u[j,i]}*fl[${j}]);
    nFr += n[${i}]*(${-u[j,i]}*fr[${j}]);

    nv += n[${i}]*${u[j,i]};
    % endfor

    nF[${j}] = nv < 0.0 ? nFl : nFr; 
% endfor
</%pyfr:macro>
