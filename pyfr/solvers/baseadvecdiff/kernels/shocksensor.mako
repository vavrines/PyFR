# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='shocksensor' ndim='1'
              du='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              artvisc='out fpdtype_t[${str(nupts)}][${str(nvars)}]'
              gradu='in fpdtype_t[${str(nupts)}][${str(nvars*ndims)}]'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'>
// Calculate grid size
fpdtype_t h = 0.0;
% for i in range(nupts):
    h += ${weights[i]}*(1.0/rcpdjac[${i}]);
% endfor
h = pow(h, ${1.0/ndims})/${order + 1};

% if sensor_type == 'rev':
    // Calculate average solution defect
    fpdtype_t int_du[${nvars}];
    fpdtype_t max_du[${nvars}];
    % for j in range(nvars):
        int_du[${j}] = 0.0;
        max_du[${j}] = 0.0;
    % endfor
    % for i,j in pyfr.ndrange(nupts, nvars):
        int_du[${j}] += ${weights[i]}*fabs(du[${i}][${j}]);
        max_du[${j}] = max(max_du[${j}], fabs(du[${i}][${j}]));
    % endfor



    % for i,j in pyfr.ndrange(nupts, nvars):
        % if vis_method == 'pointwise':
            artvisc[${i}][${j}] = ${c_mu*vis_coeffs[j]}*pow(abs(du[${i}][${j}]), ${exp_fac})*h*h*${1.0/dt_rev};
        % elif vis_method == 'max':
            artvisc[${i}][${j}] = ${c_mu*vis_coeffs[j]}*pow(max_du[${j}], ${exp_fac})*h*h*${1.0/dt_rev};
        % elif vis_method == 'mean':
            artvisc[${i}][${j}] = ${c_mu*vis_coeffs[j]}*pow(int_du[${j}], ${exp_fac})*h*h*${1.0/dt_rev};
        % elif vis_method == 'constant':
            artvisc[${i}][${j}] = ${vis_coeffs[j]*c['mu_max']};
        % endif
        artvisc[${i}][${j}] = artvisc[${i}][${j}] < ${cutoff} ? 0.0 : artvisc[${i}][${j}];
    % endfor
% elif sensor_type == 'modal':
    <% se0 = math.log10(c['s0']) %>
    // Smoothness indicator
    fpdtype_t totEn, pnEn, tmp, mu, se;
    % for var in range(nvars):
        totEn = 0.0; pnEn = 1e-15;
        % for i, deg in enumerate(ubdegs):
            tmp = ${' + '.join('{jx}*u[{j}][{svar}]'.format(j=j, jx=jx, svar=var)
                               for j, jx in enumerate(invvdm[i]) if jx != 0)};
            totEn += tmp*tmp;
            % if deg >= order:
                pnEn += tmp*tmp;
            % endif
        % endfor
        se  = ${1/math.log(10)}*log(pnEn/totEn);

        // Compute cell-wise artificial viscosity
        mu = (se < ${se0 - c['kappa']})
                     ? 0.0
                     : ${0.5*c['eps']}*h*(1.0 + sin(${0.5*math.pi/c['kappa']}*(se - ${se0})));
        mu = (se < ${se0 + c['kappa']}) ? mu : ${c['max-artvisc']};

        % for i in range(nupts):
            artvisc[${i}][${var}] = mu;
        % endfor
    % endfor
% endif 

% if scaling_type == 'ducros':
    fpdtype_t int_dscale = 0.0;
    fpdtype_t max_dscale = 0.0;
    fpdtype_t dscale;

    % if ndims == 2: 
        fpdtype_t rcprho, divrho, u_y, v_x, vort_mag_squared;
        % for i in range(nupts):
            rcprho = 1.0/u[${i}][0];
            divrho = gradu[${i}][${0*ndims + 0}] + gradu[${i}][${1*ndims + 0}];

            u_y = rcprho*(gradu[${i}][${1*ndims + 1}] - rcprho*u[${i}][1]*gradu[${i}][${1*ndims + 0}]);
            v_x = rcprho*(gradu[${i}][${0*ndims + 2}] - rcprho*u[${i}][2]*gradu[${i}][${0*ndims + 0}]);

            vort_mag_squared = pow(v_x - u_y, 2.0);

            dscale = (pow(divrho, 2.0))/(pow(divrho, 2.0) + vort_mag_squared + ${1e-6});

            int_dscale += ${weights[i]}*dscale;
            max_dscale = fmax(max_dscale, dscale);

            % if vis_method == 'pointwise':
                % for j in range(nvars):
                    artvisc[${i}][${j}] *= dscale;
                % endfor
            % endif 
        % endfor

    % elif ndims == 3:
        fpdtype_t rcprho, divrho, u_y, u_z, v_x, v_z, w_x, w_y, vort_x, vort_y, vort_z, vort_mag_squared;
        % for i in range(nupts):
            rcprho = 1/u[${i}][0];
            divrho = gradu[${i}][${0*ndims + 0}] + gradu[${i}][${1*ndims + 0}] + gradu[${i}][${2*ndims + 0}];

            u_y = rcprho*(gradu[${i}][${1*ndims + 1}] - rcprho*u[${i}][1]*gradu[${i}][${1*ndims + 0}]);
            u_z = rcprho*(gradu[${i}][${2*ndims + 1}] - rcprho*u[${i}][1]*gradu[${i}][${2*ndims + 0}]);

            v_x = rcprho*(gradu[${i}][${0*ndims + 2}] - rcprho*u[${i}][2]*gradu[${i}][${0*ndims + 0}]);
            v_z = rcprho*(gradu[${i}][${2*ndims + 2}] - rcprho*u[${i}][2]*gradu[${i}][${2*ndims + 0}]);

            w_x = rcprho*(gradu[${i}][${0*ndims + 3}] - rcprho*u[${i}][3]*gradu[${i}][${0*ndims + 0}]);
            w_y = rcprho*(gradu[${i}][${1*ndims + 3}] - rcprho*u[${i}][3]*gradu[${i}][${1*ndims + 0}]);

            vort_x = w_y - v_z;
            vort_y = u_z - w_x;
            vort_z = v_x - u_y;
            vort_mag_squared = vort_x*vort_x + vort_y*vort_y + vort_z*vort_z;

            dscale = (pow(divrho, 2.0))/(pow(divrho, 2.0) + vort_mag_squared + ${1e-6});

            int_dscale += ${weights[i]}*dscale;
            max_dscale = fmax(max_dscale, dscale);

            % if vis_method == 'pointwise':
                % for j in range(nvars):
                    artvisc[${i}][${j}] *= dscale;
                % endfor
            % endif 

        % endfor
    % endif


    % for i,j in pyfr.ndrange(nupts, nvars):
        % if vis_method == 'max':
            artvisc[${i}][${j}] *= max_dscale;
        % elif vis_method == 'mean':
            artvisc[${i}][${j}] *= int_dscale;
        % endif
    % endfor

% endif

</%pyfr:kernel>
