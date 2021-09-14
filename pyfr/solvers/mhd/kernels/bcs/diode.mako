# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mhd.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t nor = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};
    
    % for i in range(nvars):
        ur[${i}] = ul[${i}];
    % endfor

    % for i in range(ndims):
        ur[${i + 1}] -= fmin(0.0, 2*nor*nl[${i}]);
    % endfor
</%pyfr:macro>
<%pyfr:alias name='bc_rsolve_state_inv' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>


