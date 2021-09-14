# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconf_b' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='inout fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>

// Density, energy
fpdtype_t rcprho = 1.0/u[0];

% if ndims == 2:
    // Velocity and magnetic fields
    fpdtype_t v[] = {rcprho*u[1], rcprho*u[2]};
    fpdtype_t B[] = {u[3], u[4]};

    // Compute B·v
    fpdtype_t Bdotv = ${pyfr.dot('v[{i}]', 'B[{i}]', i=2)};
    fpdtype_t divB = tdivtconf[5];

    // Untransform the divergences and apply the source terms
    // Density
    tdivtconf[0] = -rcpdjac*tdivtconf[0] + ${srcex[0]};
    // Momentum
    tdivtconf[1] = -rcpdjac*(divB*B[0] + tdivtconf[1]) + ${srcex[1]};
    tdivtconf[2] = -rcpdjac*(divB*B[1] + tdivtconf[2]) + ${srcex[2]};
    // Magnetic field
    tdivtconf[3] = -rcpdjac*(divB*v[0] + tdivtconf[3]) + ${srcex[3]};
    tdivtconf[4] = -rcpdjac*(divB*v[1] + tdivtconf[4]) + ${srcex[4]};
    // DivB
    tdivtconf[5] = 0.0; 
    // Energy
    tdivtconf[6] = -rcpdjac*(divB*Bdotv + tdivtconf[6]) + ${srcex[6]}; 
    // Entropy
    tdivtconf[7] = 0.0; 

% elif ndims == 3:
    // Velocity and magnetic fields
    fpdtype_t v[] = {rcprho*u[1], rcprho*u[2], rcprho*u[3]};
    fpdtype_t B[] = {u[4], u[5], u[6]};

    // Compute B·v
    fpdtype_t Bdotv = ${pyfr.dot('v[{i}]', 'B[{i}]', i=3)};
    fpdtype_t divB = tdivtconf[7];

    // Untransform the divergences and apply the source terms
    // Density
    tdivtconf[0] = -rcpdjac*tdivtconf[0] + ${srcex[0]};
    // Momentum
    tdivtconf[1] = -rcpdjac*(divB*B[0] + tdivtconf[1]) + ${srcex[1]};
    tdivtconf[2] = -rcpdjac*(divB*B[1] + tdivtconf[2]) + ${srcex[2]};
    tdivtconf[3] = -rcpdjac*(divB*B[2] + tdivtconf[3]) + ${srcex[3]};
    // Magnetic field
    tdivtconf[4] = -rcpdjac*(divB*v[0] + tdivtconf[4]) + ${srcex[4]};
    tdivtconf[5] = -rcpdjac*(divB*v[1] + tdivtconf[5]) + ${srcex[5]};
    tdivtconf[6] = -rcpdjac*(divB*v[2] + tdivtconf[6]) + ${srcex[6]};
    // DivB
    tdivtconf[7] = 0.0; 
    // Energy
    tdivtconf[8] = -rcpdjac*(divB*Bdotv + tdivtconf[8]) + ${srcex[8]}; 
% endif

% for i in range(nvars):
    u[${i}] -= ${dt_rev}*tdivtconf[${i}];
% endfor

</%pyfr:kernel>

