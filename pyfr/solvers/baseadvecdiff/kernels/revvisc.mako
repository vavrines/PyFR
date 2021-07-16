# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='rev_viscosity_add' params='fout, rev_grads'>
% if shock_capturing == 'rev-viscosity':
% for i, j in pyfr.ndrange(ndims, nvars):
    fout[${i}][${j}] -= rev_grads[${i}][${j}];
% endfor
% endif
</%pyfr:macro>
