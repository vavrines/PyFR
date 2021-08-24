# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='ducros' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              gradu='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
              dscale='out fpdtype_t'>

% if ndims == 2:
    fpdtype_t rcprho = 1/u[0];

    fpdtype_t u_x = rcprho*(gradu[0][1] - rcprho*u[1]*gradu[0][0]);
    fpdtype_t u_y = rcprho*(gradu[1][1] - rcprho*u[1]*gradu[1][0]);
    fpdtype_t v_x = rcprho*(gradu[0][2] - rcprho*u[2]*gradu[0][0]);
    fpdtype_t v_y = rcprho*(gradu[1][2] - rcprho*u[2]*gradu[1][0]);

    fpdtype_t divu = u_x + v_y;
    fpdtype_t vort_mag_squared = pow(v_x - u_y, 2.0);
% elif ndims == 3:
    fpdtype_t rcprho = 1/u[0];

    fpdtype_t u_x = rcprho*(gradu[0][1] - rcprho*u[1]*gradu[0][0]);
    fpdtype_t u_y = rcprho*(gradu[1][1] - rcprho*u[1]*gradu[1][0]);
    fpdtype_t u_z = rcprho*(gradu[2][1] - rcprho*u[1]*gradu[2][0]);

    fpdtype_t v_x = rcprho*(gradu[0][2] - rcprho*u[2]*gradu[0][0]);
    fpdtype_t v_y = rcprho*(gradu[1][2] - rcprho*u[2]*gradu[1][0]);
    fpdtype_t v_z = rcprho*(gradu[2][2] - rcprho*u[2]*gradu[2][0]);

    fpdtype_t w_x = rcprho*(gradu[0][3] - rcprho*u[3]*gradu[0][0]);
    fpdtype_t w_y = rcprho*(gradu[1][3] - rcprho*u[3]*gradu[1][0]);
    fpdtype_t w_z = rcprho*(gradu[2][3] - rcprho*u[3]*gradu[2][0]);

    fpdtype_t divu = u_x + v_y + w_z;
    fpdtype_t vort_x = w_y - v_z;
    fpdtype_t vort_y = u_z - w_x;
    fpdtype_t vort_z = v_x - u_y;

    fpdtype_t vort_mag_squared = vort_x*vort_x + vort_y*vort_y + vort_z*vort_z;
% endif

dscale = (divu*divu)/(divu*divu + vort_mag_squared + ${1e-6});
dscale = fmax(0.0, fmin(dscale, 1.0));

</%pyfr:kernel>
