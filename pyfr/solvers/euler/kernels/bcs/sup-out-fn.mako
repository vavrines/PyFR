# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
% for i in range(nvars):
    ur[${i}] = ul[${i}];
% endfor
</%pyfr:macro>

<%pyfr:alias name='bc_rsolve_state_inv' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_rsolve_state_vis' func='bc_rsolve_state'/>