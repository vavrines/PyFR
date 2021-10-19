# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
% for i in range(ndims):
    ur[${i}] = -ul[${i}];
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
% for i in range(ndims):
    ur[${i}] = 0.0;
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, nl, grad_ul, grad_ur'>
    ${pyfr.expand('bc_common_grad_copy', 'ul', 'nl', 'grad_ul', 'grad_ur')};
</%pyfr:macro>
