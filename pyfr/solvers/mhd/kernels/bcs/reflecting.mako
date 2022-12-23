<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    ur[0] = ul[0];
    ur[${nvars-2}] = ul[${nvars-2}];
    ur[${nvars-1}] = ul[${nvars-1}];

    fpdtype_t nor = ${' + '.join(f'ul[{i + 1}]*nl[{i}]' for i in range(ndims))};
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% endfor

     nor = ${' + '.join(f'ul[{i + 1 + ndims}]*nl[{i}]' for i in range(ndims))};
% for i in range(ndims):
    ur[${i + 1 + ndims}] = ul[${i + 1 + ndims}] - 2*nor*nl[${i}];
% endfor
</%pyfr:macro>
