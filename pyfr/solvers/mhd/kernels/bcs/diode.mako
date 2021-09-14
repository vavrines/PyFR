# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mhd.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mhd.kernels.bcs.common'/>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur'>
    fpdtype_t nor = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - fmin(0.0, 2*nor*nl[${i}]);
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state_f' params='ul, gradul, artviscl, nl, magnl'>
    // Ghost state r
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve_f', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] = magnl*(ficomm[${i}]);
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state_b' params='ul, gradul, artviscl, nl, magnl'>
    // Ghost state r
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve_b', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] = magnl*(ficomm[${i}]);
% endfor
</%pyfr:macro>


<%pyfr:alias name='bc_common_flux_state_vis' func='bc_common_flux_state_f'/>
