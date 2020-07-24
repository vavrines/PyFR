# -*- coding: utf-8 -*-

<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

% if ndims == 2:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout, t, walldist'>
    fpdtype_t rho = uin[0], rhou = uin[1], rhov = uin[2], E = uin[3];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov;

    fpdtype_t rho_x = grad_uin[0][0];
    fpdtype_t rho_y = grad_uin[1][0];

    // Velocity derivatives (rho*grad[u,v])
    fpdtype_t u_x = grad_uin[0][1] - u*rho_x;
    fpdtype_t u_y = grad_uin[1][1] - u*rho_y;
    fpdtype_t v_x = grad_uin[0][2] - v*rho_x;
    fpdtype_t v_y = grad_uin[1][2] - v*rho_y;

    fpdtype_t E_x = grad_uin[0][3];
    fpdtype_t E_y = grad_uin[1][3];

% if visc_corr == 'sutherland':
    // Compute the temperature and viscosity
    fpdtype_t cpT = ${c['gamma']}*(rcprho*E - 0.5*(u*u + v*v));
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}*Trat*sqrt(Trat)
                   / (cpT + ${c['cpTs']});
% else:
    fpdtype_t mu_c = ${c['mu']};
% endif

    // Turbulence model variables and turbulent viscosity
    fpdtype_t ku = uin[4];
    fpdtype_t ku_temp = (ku > ${c['min_ku']}) ? ku : ${c['min_ku']};
    fpdtype_t wu = exp(uin[5]);

    // Turbulence model gradients 
    fpdtype_t ku_x = grad_uin[0][4];    fpdtype_t ku_y = grad_uin[1][4];
    fpdtype_t wu_x = grad_uin[0][5];    fpdtype_t wu_y = grad_uin[1][5];

	fpdtype_t sig_w2u = ${c['sig_w2']}*${c['fw']/c['fk']};
    fpdtype_t dk_dx, dw_dx, dkdw_dxi = 0;
	% for i in range(ndims): 
		dk_dx = rcprho*(grad_uin[${i}][${nvars-2}] - grad_uin[${i}][0]*uin[${nvars-2}]); 
		dw_dx = rcprho*(grad_uin[${i}][${nvars-1}] - grad_uin[${i}][0]*uin[${nvars-1}]); 
		dkdw_dxi += dk_dx*dw_dx;
	% endfor

	// Calculate damping term CDkw
	fpdtype_t CDkw = max(2*rho*${c['sig_w2']}*dkdw_dxi/wu, pow(10.0,-10));

	// Calculate blending term F1
	fpdtype_t d = walldist[0];
	fpdtype_t g1 = max(pow(ku_temp, 0.5)/(${c['betastar']}*wu*d), 500*${c['mu']}/(d*d*rho*wu));
	fpdtype_t g2 = min(g1, 4*rho*sig_w2u*ku_temp/(CDkw*d*d));
	fpdtype_t g3 = pow(g2, 4);
	fpdtype_t F1 = tanh(g3);

	// Calculate limiter term F2
	fpdtype_t g4 = max(2*pow(ku_temp, 0.5)/(${c['betastar']}*wu*d), 500*${c['mu']}/(d*d*rho*wu));
	fpdtype_t F2 = tanh(g4);

	// Calculate vorticity magnitude
	fpdtype_t vort_jk, vort_mag = 0.0, duk_dxj, duj_dxk;
	% for j,k in pyfr.ndrange(ndims,ndims):
		duk_dxj = rcprho*(grad_uin[${j}][${k+1}] - grad_uin[${j}][0]*uin[${k+1}]); // duk_dxj = 1/rho*(drhouk_dxj - drho_dxj*uk)
		duj_dxk = rcprho*(grad_uin[${k}][${j+1}] - grad_uin[${k}][0]*uin[${j+1}]); // duj_dxk = 1/rho*(drhouj_dxk - drho_dxk*uj)
		vort_jk = 0.5*(duj_dxk - duk_dxj);
		vort_mag += pow(2*vort_jk*vort_jk, 0.5);
	% endfor

    fpdtype_t sig_ku = ${c['fw']/c['fk']}*(F1*${c['sig_k1']} + (1-F1)*${c['sig_k2']});
    fpdtype_t sig_wu = ${c['fw']/c['fk']}*(F1*${c['sig_w1']} + (1-F1)*${c['sig_w2']});

    fpdtype_t mu_t = rho*${c['a1']}*ku/max(${c['a1']}*wu, vort_mag*F2);

	mu_t = (mu_t < 0.0) ? 0.0 : mu_t;
	mu_t = ${c['tmswitch']}*(1.0 - exp(-${c['tdvc']}*(t - ${c['tmstarttime']})))*mu_t;
	mu_t = (mu_t > ${c['mu']}*${c['max_mutrat']}) ? ${c['mu']}*${c['max_mutrat']} : mu_t;


    // Compute temperature derivatives (c_v*dT/d[x,y])
    fpdtype_t T_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x));
    fpdtype_t T_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y));

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*(mu_c + mu_t)*rcprho*(u_x - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_yy = -2*(mu_c + mu_t)*rcprho*(v_y - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_xy = -(mu_c + mu_t)*rcprho*(v_x + u_y);

    fout[0][1] += t_xx;     fout[1][1] += t_xy;
    fout[0][2] += t_xy;     fout[1][2] += t_yy;

    fout[0][3] += u*t_xx + v*t_xy + -(mu_c*${c['gamma']/c['Pr']} + mu_t*${c['gamma']/c['Pr_t']})*T_x;
    fout[1][3] += u*t_xy + v*t_yy + -(mu_c*${c['gamma']/c['Pr']} + mu_t*${c['gamma']/c['Pr_t']})*T_y;

    fout[0][4] += -rcprho*(mu_c + sig_ku*mu_t)*ku_x;     fout[1][4] += -rcprho*(mu_c + sig_ku*mu_t)*ku_y; 
    fout[0][5] += -rcprho*(mu_c + sig_wu*mu_t)*wu_x;     fout[1][5] += -rcprho*(mu_c + sig_wu*mu_t)*wu_y; 



</%pyfr:macro>
% elif ndims == 3:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout, t, walldist'>
    fpdtype_t rho  = uin[0];
    fpdtype_t rhou = uin[1], rhov = uin[2], rhow = uin[3];
    fpdtype_t E    = uin[4];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov, w = rcprho*rhow;

    fpdtype_t rho_x = grad_uin[0][0];
    fpdtype_t rho_y = grad_uin[1][0];
    fpdtype_t rho_z = grad_uin[2][0];

    // Velocity derivatives (rho*grad[u,v,w])
    fpdtype_t u_x = grad_uin[0][1] - u*rho_x;
    fpdtype_t u_y = grad_uin[1][1] - u*rho_y;
    fpdtype_t u_z = grad_uin[2][1] - u*rho_z;
    fpdtype_t v_x = grad_uin[0][2] - v*rho_x;
    fpdtype_t v_y = grad_uin[1][2] - v*rho_y;
    fpdtype_t v_z = grad_uin[2][2] - v*rho_z;
    fpdtype_t w_x = grad_uin[0][3] - w*rho_x;
    fpdtype_t w_y = grad_uin[1][3] - w*rho_y;
    fpdtype_t w_z = grad_uin[2][3] - w*rho_z;

    fpdtype_t E_x = grad_uin[0][4];
    fpdtype_t E_y = grad_uin[1][4];
    fpdtype_t E_z = grad_uin[2][4];

% if visc_corr == 'sutherland':
    // Compute the temperature and viscosity
    fpdtype_t cpT = ${c['gamma']}*(rcprho*E - 0.5*(u*u + v*v + w*w));
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}*Trat*sqrt(Trat)
                   / (cpT + ${c['cpTs']});
% else:
    fpdtype_t mu_c = ${c['mu']};
% endif

    fpdtype_t ku = uin[5];
    fpdtype_t ku_temp = (ku > ${c['min_ku']}) ? ku : ${c['min_ku']};
    fpdtype_t wu = exp(uin[5]);

    // Turbulence model gradients 
    fpdtype_t ku_x = grad_uin[0][5];    fpdtype_t ku_y = grad_uin[1][5];    fpdtype_t ku_z = grad_uin[2][5];
    fpdtype_t wu_x = grad_uin[0][6];    fpdtype_t wu_y = grad_uin[1][6];    fpdtype_t wu_z = grad_uin[2][6];

	fpdtype_t sig_w2u = ${c['sig_w2']}*${c['fw']/c['fk']};
    fpdtype_t dk_dx, dw_dx, dkdw_dxi = 0;
	% for i in range(ndims): 
		dk_dx = rcprho*(grad_uin[${i}][${nvars-2}] - grad_uin[${i}][0]*uin[${nvars-2}]); 
		dw_dx = rcprho*(grad_uin[${i}][${nvars-1}] - grad_uin[${i}][0]*uin[${nvars-1}]); 
		dkdw_dxi += dk_dx*dw_dx;
	% endfor

	// Calculate damping term CDkw
	fpdtype_t CDkw = max(2*rho*${c['sig_w2']}*dkdw_dxi/wu, pow(10.0,-10));

	// Calculate blending term F1
	fpdtype_t d = walldist[0];
	fpdtype_t g1 = max(pow(ku_temp, 0.5)/(${c['betastar']}*wu*d), 500*${c['mu']}/(d*d*rho*wu));
	fpdtype_t g2 = min(g1, 4*rho*sig_w2u*ku_temp/(CDkw*d*d));
	fpdtype_t g3 = pow(g2, 4);
	fpdtype_t F1 = tanh(g3);

	// Calculate limiter term F2
	fpdtype_t g4 = max(2*pow(ku_temp, 0.5)/(${c['betastar']}*wu*d), 500*${c['mu']}/(d*d*rho*wu));
	fpdtype_t F2 = tanh(g4);

	// Calculate vorticity magnitude
	fpdtype_t vort_jk, vort_mag = 0.0, duk_dxj, duj_dxk;
	% for j,k in pyfr.ndrange(ndims,ndims):
		duk_dxj = rcprho*(grad_uin[${j}][${k+1}] - grad_uin[${j}][0]*uin[${k+1}]); // duk_dxj = 1/rho*(drhouk_dxj - drho_dxj*uk)
		duj_dxk = rcprho*(grad_uin[${k}][${j+1}] - grad_uin[${k}][0]*uin[${j+1}]); // duj_dxk = 1/rho*(drhouj_dxk - drho_dxk*uj)
		vort_jk = 0.5*(duj_dxk - duk_dxj);
		vort_mag += pow(2*vort_jk*vort_jk, 0.5);
	% endfor

    fpdtype_t sig_ku = ${c['fw']/c['fk']}*(F1*${c['sig_k1']} + (1-F1)*${c['sig_k2']});
    fpdtype_t sig_wu = ${c['fw']/c['fk']}*(F1*${c['sig_w1']} + (1-F1)*${c['sig_w2']});

    fpdtype_t mu_t = rho*${c['a1']}*ku/max(${c['a1']}*wu, vort_mag*F2);

	mu_t = (mu_t < 0.0) ? 0.0 : mu_t;
	mu_t = ${c['tmswitch']}*(1.0 - exp(-${c['tdvc']}*(t - ${c['tmstarttime']})))*mu_t;
	mu_t = (mu_t > ${c['mu']}*${c['max_mutrat']}) ? ${c['mu']}*${c['max_mutrat']} : mu_t;

    // Compute temperature derivatives (c_v*dT/d[x,y,z])
    fpdtype_t T_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x + w*w_x));
    fpdtype_t T_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y + w*w_y));
    fpdtype_t T_z = rcprho*(E_z - (rcprho*rho_z*E + u*u_z + v*v_z + w*w_z));

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*(mu_c + mu_t)*rcprho*(u_x - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_yy = -2*(mu_c + mu_t)*rcprho*(v_y - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_zz = -2*(mu_c + mu_t)*rcprho*(w_z - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_xy = -(mu_c + mu_t)*rcprho*(v_x + u_y);
    fpdtype_t t_xz = -(mu_c + mu_t)*rcprho*(u_z + w_x);
    fpdtype_t t_yz = -(mu_c + mu_t)*rcprho*(w_y + v_z);

    fout[0][1] += t_xx;     fout[1][1] += t_xy;     fout[2][1] += t_xz;
    fout[0][2] += t_xy;     fout[1][2] += t_yy;     fout[2][2] += t_yz;
    fout[0][3] += t_xz;     fout[1][3] += t_yz;     fout[2][3] += t_zz;

    fout[0][4] += u*t_xx + v*t_xy + w*t_xz + -(mu_c*${c['gamma']/c['Pr']} + mu_t*${c['gamma']/c['Pr_t']})*T_x;
    fout[1][4] += u*t_xy + v*t_yy + w*t_yz + -(mu_c*${c['gamma']/c['Pr']} + mu_t*${c['gamma']/c['Pr_t']})*T_y;
    fout[2][4] += u*t_xz + v*t_yz + w*t_zz + -(mu_c*${c['gamma']/c['Pr']} + mu_t*${c['gamma']/c['Pr_t']})*T_z;

    fout[0][5] += -rcprho*(mu_c + sig_ku*mu_t)*ku_x;     fout[1][5] += -rcprho*(mu_c + sig_ku*mu_t)*ku_y;     fout[2][5] += -rcprho*(mu_c + sig_ku*mu_t)*ku_z; 
    fout[0][6] += -rcprho*(mu_c + sig_wu*mu_t)*wu_x;     fout[1][6] += -rcprho*(mu_c + sig_wu*mu_t)*wu_y;     fout[2][6] += -rcprho*(mu_c + sig_wu*mu_t)*wu_z; 

</%pyfr:macro>
% endif
