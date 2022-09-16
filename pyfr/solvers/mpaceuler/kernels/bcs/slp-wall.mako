# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, nl, ur, qr'>
    fpdtype_t nor = ${' + '.join('nl[{0}]*ul[{1}]'.format(i, i + 1)
                                 for i in range(ndims))};
    ur[0] = ul[0];
    ur[${ndims + 1}] = ul[${ndims + 1}];
    qr[0] = ql[0];
    
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% endfor
</%pyfr:macro>
