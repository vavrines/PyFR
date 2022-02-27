<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.mhd.kernels.rsolvers.${rsolver}'/>

<% tau = c['ldg-tau'] %>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, artviscl, nl, magnl'>
    // Inviscid (Riemann solve) state
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}], fvcomm;
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

    % for i in range(nvars):
        ul[${i}] = magnl*ficomm[${i}];
    % endfor
</%pyfr:macro>

