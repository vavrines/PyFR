# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t nor = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};
// Mass
% for i in range(nvars):
    ur[${i}] = ul[${i}];
% endfor
// Momentum
% for i in range(ndims):
    ur[${i + nspec}] = ul[${i + nspec}] - 2*nor*nl[${i}];
% endfor
// Energy + fractions
% for i in range(nspec):
    ur[${i + nspec + ndims}] = ul[${i + nspec + ndims}];
% endfor
</%pyfr:macro>
