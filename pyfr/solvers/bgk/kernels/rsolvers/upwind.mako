# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='rsolve' params='fl, fr, n, nF, u'>
    fpdtype_t nv;

for (int j = 0; j < ${nvars}; j++)
{
    nF[j] = nv = 0.0;
    % for i in range(ndims):
    nv += n[${i}]*u[j][${i}];
    % endfor

    if (nv > 0.0) {
        % for i in range(ndims):
        nF[j] += n[${i}]*(u[j][${i}]*fl[j]);
        % endfor
    }
    else {
        % for i in range(ndims):
        nF[j] += n[${i}]*(u[j][${i}]*fr[j]);
        % endfor
    }
}
</%pyfr:macro>
