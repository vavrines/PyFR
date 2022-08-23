# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr, u, M' externs='ploc, t'>
% for i in range(nvars):
    fr[${i}] = fl[${i}];
% endfor
</%pyfr:macro>
