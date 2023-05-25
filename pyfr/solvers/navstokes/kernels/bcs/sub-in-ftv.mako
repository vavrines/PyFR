<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t pl = ${c['gamma'] - 1.0}*(ul[${nvars - 1}]
                 - (0.5/ul[0])*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))});
    ur[0] = pl/${c['theta']};
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = ur[0] * (${c[v]});
% endfor
    ur[${nvars - 1}] = ${1.0/(c['gamma'] - 1.0)}*pl + 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
