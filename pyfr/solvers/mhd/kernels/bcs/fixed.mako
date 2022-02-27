# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mhd.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mhd.kernels.bcs.common'/>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur'>

    % if ndims == 2:
        ur[0] = ${c['rho']};
        ur[1] = (${c['rho']})*(${c['u']});
        ur[2] = (${c['rho']})*(${c['v']});
        ur[3] = ${c['Bx']};
        ur[4] = ${c['By']};
        ur[5] = 0.0;
        ur[6] = ${c['p']}/(${c['gamma']-1.}) + (0.5/(${c['rho']}))*(ur[1]*ur[1] + ur[2]*ur[2]) 
                                             + 0.5*(ur[3]*ur[3] + ur[4]*ur[4]);
    % elif ndims == 3:
        ur[0] = ${c['rho']};
        ur[1] = (${c['rho']})*(${c['u']});
        ur[2] = (${c['rho']})*(${c['v']});
        ur[3] = (${c['rho']})*(${c['w']});
        ur[4] = ${c['Bx']};
        ur[5] = ${c['By']};
        ur[6] = ${c['Bz']};
        ur[7] = 0.0;
        ur[8] = ${c['p']}/(${c['gamma']-1.}) + (0.5/(${c['rho']}))*(ur[1]*ur[1] + ur[2]*ur[2] + ur[3]*ur[3]) 
                                             + 0.5*(ur[4]*ur[4] + ur[5]*ur[5] + ur[6]*ur[6]);
    % endif

</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, artviscl, nl, magnl'>
    // Ghost state r
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] = magnl*(ficomm[${i}]);
% endfor
</%pyfr:macro>
