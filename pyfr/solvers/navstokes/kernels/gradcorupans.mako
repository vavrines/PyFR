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
              walldist='in fpdtype_t[${str(ndims)}]'>


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
	prod += duj_dxk*Tjk; // CHECK THESIS?
% endfor


// Multiplying by sig_w2u instead of dividing so fw/fk not fk/fw
fpdtype_t sig_w2u = ${c['sig_w2']}*${c['fw']/c['fk']};

// Production limiter (Menter, F. R., AIAA Paper 93-2906, July 1993)
//prod = min(prod, 20*${c['betastar']}*rho*wu*ku_temp);

// Convert to unresolved production
fpdtype_t prod_u = ${c['fk']}*prod + ${c['betastar']}*ku_temp*wu*(1.0 - 1.0/${c['fw']});

// Calculate damping term CDkw
fpdtype_t CDkw = max(2*rho*sig_w2u*dkdw_dxi/wu, pow(10.0,-10));


// Calculate blending term F1
fpdtype_t d = walldist[0];
fpdtype_t g1 = max(pow(ku_temp, 0.5)/(${c['betastar']}*wu*d), 500*${c['mu']}/(d*d*rho*wu));
fpdtype_t g2 = min(g1, 4*rho*sig_w2u*ku_temp/(CDkw*d*d));
fpdtype_t g3 = pow(g2, 4);
fpdtype_t F1 = tanh(g3);

// Calculate blended constants
fpdtype_t alpha = F1*${c['alpha1']} + (1 - F1)*${c['alpha2']};
fpdtype_t beta  = F1*${c['beta1']}  + (1 - F1)*${c['beta2']};
fpdtype_t betaprime = alpha*${c['betastar']} - alpha*${c['betastar']}/${c['fw']} + beta/${c['fw']};



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
