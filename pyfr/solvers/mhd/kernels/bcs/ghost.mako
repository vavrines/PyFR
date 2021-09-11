<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.mhd.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mhd.kernels.flux'/>

<% tau = c['ldg-tau'] %>

<%pyfr:macro name='bc_common_flux_state_f' params='ul, gradul, artviscl, nl, magnl'>
    // Inviscid (Riemann solve) state
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state_inv', 'ul', 'nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}], fvcomm;
    ${pyfr.expand('rsolve_f', 'ul', 'ur', 'nl', 'ficomm')};

    % for i in range(nvars):
        ul[${i}] = magnl*ficomm[${i}];
    % endfor
</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state_b' params='ul, gradul, artviscl, nl, magnl'>
    // Inviscid (Riemann solve) state
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state_inv', 'ul', 'nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}], fvcomm;
    ${pyfr.expand('rsolve_b', 'ul', 'ur', 'nl', 'ficomm')};

    % for i in range(nvars):
        ul[${i}] = magnl*ficomm[${i}];
    % endfor
</%pyfr:macro>
