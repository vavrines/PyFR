# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='artificial_viscosity_add' params='grad_uin, u, fout, artvisc'>
% if shock_capturing == 'artificial-viscosity':
    fpdtype_t scale = 1.0;
    % if scaling_type == 'ducros':
        % if ndims == 2:
            fpdtype_t rcprho = 1/u[0];
            fpdtype_t divrho = grad_uin[0][0] + grad_uin[1][0];

            fpdtype_t u_y = rcprho*(grad_uin[1][1] - rcprho*u[1]*grad_uin[1][0]);
            fpdtype_t v_x = rcprho*(grad_uin[0][2] - rcprho*u[2]*grad_uin[0][0]);

            fpdtype_t vort_mag_squared = pow(v_x - u_y, 2.0);
        % elif ndims == 3:
            fpdtype_t rcprho = 1/u[0];
            fpdtype_t divrho = grad_uin[0][0] + grad_uin[1][0] + grad_uin[2][0];

            fpdtype_t u_y = rcprho*(grad_uin[1][1] - rcprho*u[1]*grad_uin[1][0]);
            fpdtype_t u_z = rcprho*(grad_uin[2][1] - rcprho*u[1]*grad_uin[2][0]);

            fpdtype_t v_x = rcprho*(grad_uin[0][2] - rcprho*u[2]*grad_uin[0][0]);
            fpdtype_t v_z = rcprho*(grad_uin[2][2] - rcprho*u[2]*grad_uin[2][0]);

            fpdtype_t w_x = rcprho*(grad_uin[0][3] - rcprho*u[3]*grad_uin[0][0]);
            fpdtype_t w_y = rcprho*(grad_uin[1][3] - rcprho*u[3]*grad_uin[1][0]);

            fpdtype_t vort_x = w_y - v_z;
            fpdtype_t vort_y = u_z - w_x;
            fpdtype_t vort_z = v_x - u_y;

            fpdtype_t vort_mag_squared = vort_x*vort_x + vort_y*vort_y + vort_z*vort_z;
        % endif

        scale = (divrho*divrho)/(divrho*divrho + vort_mag_squared + ${1e-6});
    % endif 

    % for i, j in pyfr.ndrange(ndims, nvars):
        fout[${i}][${j}] -= fmin(scale*artvisc[${j}], ${mu_max})*grad_uin[${i}][${j}];
    % endfor
% endif

</%pyfr:macro>
