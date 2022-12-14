<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t nor = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};
    
    % for i in range(nvars):
    ur[${i}] = ul[${i}];
    % endfor

    if (nor < 0.0) {
        % for i in range(ndims):
        ur[${i + 1}] -= 2*nor*nl[${i}];
        % endfor
    }
</%pyfr:macro>
