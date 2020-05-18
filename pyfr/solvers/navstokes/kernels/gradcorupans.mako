# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='gradcorupans' ndim='2'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              rcpdjac='in fpdtype_t'
              gradu='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'
              u='in fpdtype_t[${str(nvars)}]'
              ku_src='inout fpdtype_t'
              wu_src='inout fpdtype_t'
              t = 'scalar fpdtype_t'
              ploc = 'in fpdtype_t[${str(ndims)}]'
              F1='inout fpdtype_t'
              fk='in fpdtype_t'>

fpdtype_t tmpgradu[${ndims}];

// Get density gradients first
% for i in range(ndims):
    tmpgradu[${i}] = gradu[${i}][${0}];
% endfor	
% for i in range(ndims):
    gradu[${i}][${0}] = rcpdjac*(${' + '.join('smats[{k}][{i}]*tmpgradu[{k}]'
                                              .format(i=i, k=k)
                                              for k in range(ndims))});
% endfor


% for j in range(ndims):
	% for i in range(ndims):
	    tmpgradu[${i}] = gradu[${i}][${j+1}];
	% endfor

	% for i in range(ndims):
	    gradu[${i}][${j+1}] = rcpdjac*(${' + '.join('smats[{k}][{i}]*tmpgradu[{k}]'
	                                              .format(i=i, k=k)
	                                              for k in range(ndims))});
	% endfor

% endfor

// Setup variables
fpdtype_t prod = 0.0;
fpdtype_t rho = u[0];
fpdtype_t rcprho = 1/u[0];
fpdtype_t duk_dxj, duj_dxk;

fpdtype_t ku = u[${nvars-2}];
fpdtype_t wu = exp(u[${nvars-1}]);

// Mu_t must be positive
fpdtype_t mu_t = (rho*ku/wu < 0.0) ? 0.0 : rho*ku/wu;

// Adding viscosity rampup and max mu_t ratio limiter
mu_t = ${c['tmswitch']}*(1.0 - exp(-${c['tdvc']}*(t - ${c['tmstarttime']})))*mu_t;
mu_t = (mu_t > ${c['mu']}*${c['max_mutrat']}) ? ${c['mu']}*${c['max_mutrat']} : mu_t;


// Get production term
fpdtype_t Sjk = 0.0;
fpdtype_t Tjk = 0.0;
fpdtype_t dui_dxi = 0.0;
fpdtype_t dk_dx, dw_dx;
fpdtype_t dkdw_dxi = 0.0;

fpdtype_t ku_temp = (ku < ${c['min_ku']}) ? ${c['min_ku']} : ku;

% for i in range(ndims):
	dui_dxi += rcprho*(gradu[${i}][${i+1}] - gradu[${i}][0]*u[${i+1}]); 
	dk_dx = rcprho*(gradu[${i}][${nvars-2}] - gradu[${i}][0]*u[${nvars-2}]); 
	dw_dx = rcprho*(gradu[${i}][${nvars-1}] - gradu[${i}][0]*u[${nvars-1}]); 
	dkdw_dxi += dk_dx*dw_dx;
% endfor

% for j,k in pyfr.ndrange(ndims,ndims):
	duk_dxj = rcprho*(gradu[${j}][${k+1}] - gradu[${j}][0]*u[${k+1}]); // duk_dxj = 1/rho*(drhouk_dxj - drho_dxj*uk)
	duj_dxk = rcprho*(gradu[${k}][${j+1}] - gradu[${k}][0]*u[${j+1}]); // duj_dxk = 1/rho*(drhouj_dxk - drho_dxk*uj)

	Sjk = 0.5*(duk_dxj + duj_dxk);
	% if (j == k):
		Tjk = mu_t*(2*Sjk - ${2.0/3.0}*dui_dxi) - ${2.0/3.0}*rho*ku_temp;
	% else:
		Tjk = mu_t*(2*Sjk);
	% endif
	prod += duj_dxk*Tjk; 
% endfor


fpdtype_t fk_temp;
% if c['BLS'] > 0.5:
	fk_temp = (1 - F1)*fk; // Boundary layer shielding
% else:
	fk_temp = fk;
% endif

fk_temp = min(${c['max_fk']}, max(${c['min_fk']}, fk_temp));
fpdtype_t fw = 1.0/fk_temp; // Assume fw = 1/fk

fpdtype_t sig_w2u = ${c['sig_w2']}*fw/fk_temp;

// Production limiter (Menter, F. R., AIAA Paper 93-2906, July 1993)
//prod = min(prod, 20*${c['betastar']}*rho*wu*ku_temp);

// Convert to unresolved production
fpdtype_t prod_u = fk_temp*prod + ${c['betastar']}*ku_temp*wu*(1.0 - 1.0/fw);

// Calculate damping term CDkw
fpdtype_t CDkw = max(2*rho*${c['sig_w2']}*dkdw_dxi/wu, pow(10.0,-10));

// d = sqrt(x**2 + y**2) - 0.5 for cylinder of diameter 1
fpdtype_t d;
% if geo == 'cylinder':
	d = pow(pow(ploc[0], 2) + pow(ploc[1], 2), 0.5) - 0.5; // Cylinder
% elif geo == 'tandsphere':
	d = min(pow(pow(ploc[0], 2) + pow(ploc[1], 2), 0.5) - 0.5, pow(pow(ploc[0]-10, 2) + pow(ploc[1], 2), 0.5) - 0.5); // Tandem spheres
% elif geo == 'bfstep':
	if (ploc[0] > 0.0 && ploc[0] < 1.0 && ploc[1] > 1.0){
		d = pow(pow(ploc[0], 2) + pow(ploc[1]-1, 2), 0.5);
	}
	else {
		d = (ploc[0] <= 0.0) ? ploc[1] - 1.0 : ploc[1];
	}
% elif geo == 'cube':
	fpdtype_t d1 = max(0.0, abs(ploc[0]) - 0.5);
	fpdtype_t d2 = max(0.0, abs(ploc[1]) - 1.0);
	fpdtype_t d3 = max(0.0, abs(ploc[2]) - 0.5);
	d = pow(pow(d1,2) + pow(d2,2) + pow(d3,2), 0.5);
	d = min(d, ploc[1]);
% elif geo == 'TGV':
	d = 1000000;
% endif

// Calculate blending term F1
fpdtype_t g1 = max(pow(ku_temp, 0.5)/(${c['betastar']}*wu*d), 500*${c['mu']}/(d*d*rho*wu));
fpdtype_t g2 = min(g1, 4*rho*sig_w2u*ku_temp/(CDkw*d*d));
fpdtype_t g3 = pow(g2, 4);
F1 = tanh(g3);

// Calculate blended constants
fpdtype_t alpha = F1*${c['alpha1']} + (1 - F1)*${c['alpha2']};
fpdtype_t beta  = F1*${c['beta1']}  + (1 - F1)*${c['beta2']};
fpdtype_t betaprime = alpha*${c['betastar']} - alpha*${c['betastar']}/fw + beta/fw;



// Calculate ku and wu source terms

ku_src = ${c['tmswitch']}*(prod_u - ${c['betastar']}*rho*ku_temp*wu);
wu_src = ${c['tmswitch']}*(alpha*prod_u*wu/ku_temp - rho*betaprime*wu*wu) + 2*(1-F1)*(rho*sig_w2u/wu)*dkdw_dxi;

ku_src = (ku < ${c['min_ku']}) ? ${c['ku_limiter']} : ku_src;


// Get gradients for energy and turbulence model variables
% for j in range(ndims+1, nvars):
	% for i in range(ndims):
	    tmpgradu[${i}] = gradu[${i}][${j}];
	% endfor
	% for i in range(ndims):
	    gradu[${i}][${j}] = rcpdjac*(${' + '.join('smats[{k}][{i}]*tmpgradu[{k}]'
	                                              .format(i=i, k=k)
	                                              for k in range(ndims))});
		
	% endfor
% endfor


</%pyfr:kernel>
