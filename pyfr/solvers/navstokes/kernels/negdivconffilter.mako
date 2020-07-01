# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconffilter' ndim='1'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ploc='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t'
              shockcell='in fpdtype_t'>

fpdtype_t solpred[${nupts}][${nvars}], magratios[${nupts-1}];
fpdtype_t solpredmodes[${nupts}][${nvars}], filteredsolmodes[${nupts}][${nvars}], filtercoeffs[${nupts}];
fpdtype_t tmp;

// Calculate transformed divF
% for i in range(nupts):
	% for j, ex in enumerate(srcex):
	    tdivtconf[${i}][${j}] = -rcpdjac*tdivtconf[${i}][${j}] + ${ex};
	% endfor
% endfor

// For cells with shocks, apply filtering operation
// Solution filtering is performed by modifying the divergence of the flux such that the solution at the next time step 
// 		is the filtered solution

// To avoid compiler warnings, only declare variables if needed for filter
% if filtermethod == 'linear':
	fpdtype_t maxslope = 0.0, slope;
% elif filtermethod == 'exponential':
	fpdtype_t maxalpha = 0.0, alpha;
% endif

if (shockcell == 1){

	// Loop to individually filter each variable
	% for svar in range(nvars):
		// Calculate solpred = forward Euler prediction of u[shockvar] at the next time step
		% for i in range(nupts):
			solpred[${i}][${svar}] = u[${i}][${svar}] + ${dt}*tdivtconf[${i}][${svar}];
		% endfor

		// Calculate solpredmodes = modal coefficients of solpred
		% for i in range(nupts):
			tmp = ${' + '.join('{jx}*solpred[{j}][{k}]'.format(j=j, jx=jx, k=svar) for j, jx in enumerate(invvdm[i]) if jx != 0)};
			solpredmodes[${i}][${svar}] = tmp;
		% endfor

		// Calculate ratio of solution magnitudes to target ratio
		% for i in range(nupts-1):
			magratios[${i}] = pow(solpredmodes[${i+1}][${svar}]/solpredmodes[0][${svar}], 2.0)/${targetratios[i]};
		% endfor

		// Filter using pointwise, linear, or exponential filter
		// Pointwise: filters modes above decay rate to the decay rate
		// Linear: filters modes with a linear slope determined by the mode magnitude most exceeding decay rate
		// Exponential: filters modes with a exponential decay determined by the mode magnitude most exceeding decay rate
		//				Exponential filter coefficient B in exp(-a*k*B) is a parameter (suitable range 1 to 4)

		filtercoeffs[0] = 1.0; // Leave mean mode unfiltered
		% if filtermethod == 'pointwise':
			// Compute filtercoeffs = pointwise filter coefficients
			% for i in range(nupts-1):
				if (magratios[${i}] > 1.0){
					filtercoeffs[${i+1}] = fmax(0.0, 1.0/magratios[${i}]);
				}
				else {
					filtercoeffs[${i+1}] = 1.0;
				}
			% endfor

		% elif filtermethod == 'linear':
			maxslope = 0.0;
			% for i in range(nupts-1):			
				if (magratios[${i}] > 1.0){
					// Find slope needed such that this mode fits within the decay rate
					slope = (magratios[${i}] - 1.0)/${ubdegs[i+1]};
					// Find largest slope needed
					if (slope > maxslope){ 
						maxslope = slope;
					}
				}
			% endfor

			% for i in range(nupts-1):
				filtercoeffs[${i+1}] = fmax(0.0, 1.0/(1.0 + maxslope*${ubdegs[i+1]})); // Set filter to 0 if negative
			% endfor

		% elif filtermethod == 'exponential':
			maxalpha = 0.0;
			% for i in range(nupts-1):			
				if (magratios[${i}] > 1.0){
					// Find exponent factor alpha such that this mode fits within the decay rate
					alpha = log(1.0/magratios[${i}])/(-pow(${(i+1.)/max(ubdegs)}, ${kexp}));
					// Find largest factor needed
					if (alpha > maxalpha){ 
						maxalpha = alpha;
					}
				}
			% endfor

			% for i in range(nupts-1):
				filtercoeffs[${i+1}] = exp(-maxalpha*pow(${(i+1.)/max(ubdegs)}, ${kexp}));
			% endfor
		% endif

		// Compute filteredsolmodes = filtered predicted solution modes 
		% for i in range(nupts):
			filteredsolmodes[${i}][${svar}] = filtercoeffs[${i}]*solpredmodes[${i}][${svar}];
		% endfor


		// Compute filteredsol = nodal form of filtered solution modes 
		// Compute divF such that the solution at the next time step is the filtered solution
		% for i in range(nupts):
			tmp = ${' + '.join('{jx}*filteredsolmodes[{j}][{k}]'.format(j=j, jx=jx, k=svar) for j, jx in enumerate(vdm[i]) if jx != 0)};
			tdivtconf[${i}][${svar}] = (tmp - u[${i}][${svar}])/${dt};
		% endfor
	% endfor
}




</%pyfr:kernel>
