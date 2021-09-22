# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mhd.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    % for i in range(nvars):
        ur[${i}] = ul[${i}];
    % endfor
    fpdtype_t norv = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};
    fpdtype_t norb = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1 + ndims)
                                 for i in range(ndims))};
    % for i in range(ndims):
        ur[${i + 1}] = ul[${i + 1}] - 2*norv*nl[${i}];
        ur[${i + 1 + ndims}] = ul[${i + 1 + ndims}] - 2*norb*nl[${i}];
    % endfor

</%pyfr:macro>
<%pyfr:alias name='bc_rsolve_state_inv' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
