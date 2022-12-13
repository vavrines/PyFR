<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mhd.kernels.rsolvers.${rsolver}'/>

<% tol = 1e-6 %>

// Performs rotation from unit vector n to unit vector m for vector v
<%pyfr:macro name='rotate' params='v, n, m'>
    // Rotation matrix 
    // [ a   b
    //   c   d ]
    fpdtype_t a, b, c, d;
    a = n[0]*m[0] + n[1]*m[1];
    b = m[0]*n[1] - n[0]*m[1];
    c = n[0]*m[1] - m[0]*n[1];
    d = n[0]*m[0] + n[1]*m[1];

    // Rotate only u/v and Bx/By components
    fpdtype_t v1 = v[1];
    fpdtype_t v2 = v[2];
    v[1] = a*v1 + b*v2;
    v[2] = c*v1 + d*v2;

    fpdtype_t b1 = v[${ndims+1}];
    fpdtype_t b2 = v[${ndims+2}];
    v[${ndims+1}] = a*b1 + b*b2;
    v[${ndims+2}] = c*b1 + d*b2;

</%pyfr:macro>

<%pyfr:kernel name='intcflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='inout view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              nr='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_n = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});

    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_n)*nl[{i}]', i=ndims)};
    fpdtype_t negnorm_nr[] = ${pyfr.array('-(1 / mag_n)*nr[{i}]', i=ndims)};

    // If rotationally periodic, perform transformations 
    if ( (abs(norm_nl[0] - negnorm_nr[0]) > ${tol}) || (abs(norm_nl[1] - negnorm_nr[1]) > ${tol}) ) {
        ${pyfr.expand('rotate', 'ur', 'negnorm_nr', 'norm_nl')};
        // Perform the Riemann solve
        fpdtype_t fn[${nvars}];

        ${pyfr.expand('rsolve', 'ul', 'ur', 'norm_nl', 'fn')};
        % for i in range(nvars):
        ul[${i}] =  mag_n*(fn[${i}]);
        ur[${i}] = -mag_n*(fn[${i}]);
        % endfor

        // Transform RHS flux
        ${pyfr.expand('rotate', 'ur', 'norm_nl', 'negnorm_nr')};
    }
    else {
        // Perform the Riemann solve
        fpdtype_t fn[${nvars}];
        ${pyfr.expand('rsolve', 'ul', 'ur', 'norm_nl', 'fn')};

        // Scale and write out the common normal fluxes
        % for i in range(nvars):
        ul[${i}] =  mag_n*fn[${i}];
        ur[${i}] = -mag_n*fn[${i}];
        % endfor
    }
</%pyfr:kernel>
