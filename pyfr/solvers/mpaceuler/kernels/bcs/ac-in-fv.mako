# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, nl, ur, qr' externs='ploc, t'>
    ur[0] = ul[0];
    ur[${ndims + 1}] = ul[${ndims + 1}];
    qr[0] = ql[0];

% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = ql[0]*${c[v]};
% endfor
</%pyfr:macro>
