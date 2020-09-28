# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t nor = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};
    ur[0] = ul[0];
    fpdtype_t tmp;
% for i in range(ndims):
	tmp = ul[${i + 1}] - 2*nor*nl[${i}];
    ul[${i + 1}] = tmp;
    ur[${i + 1}] = tmp;
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>
