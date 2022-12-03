<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

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

    // Offset by 1 (rotate only u/v velocity components)
    fpdtype_t v1 = v[1];
    fpdtype_t v2 = v[2];

    v[1] = a*v1 + b*v2;
    v[2] = c*v1 + d*v2;
</%pyfr:macro>

<%pyfr:kernel name='pintconu' ndim='1'
              ulin='in view fpdtype_t[${str(nvars)}]'
              urin='in view fpdtype_t[${str(nvars)}]'
              ulout='out view fpdtype_t[${str(nvars)}]'
              urout='out view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              nr='in fpdtype_t[${str(ndims)}]'>

    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    fpdtype_t mag_nr = sqrt(${pyfr.dot('nr[{i}]', i=ndims)});
    fpdtype_t negnorm_nr[] = ${pyfr.array('-(1 / mag_nr)*nr[{i}]', i=ndims)};

    fpdtype_t tmpu[${nvars}];

    // ----- Assumes c['ldg-beta'] == 0.5 -----
    % for i in range(nvars):
    tmpu[${i}] = urin[${i}];
    % endfor
    ${pyfr.expand('rotate', 'tmpu', 'negnorm_nr', 'norm_nl')};
    % for i in range(nvars):
    ulout[${i}] = urin[${i}];
    % endfor
    // ----------------------------------------

</%pyfr:kernel>
