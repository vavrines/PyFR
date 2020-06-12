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

fpdtype_t solpred[${nupts}][${nvars}],      magratios[${nupts-1}];
fpdtype_t solpredmodes[${nupts}][${nvars}], filteredsolmodes[${nupts}][${nvars}], filtercoeffs[${nupts}];
fpdtype_t tmp, targetratio;

// Calculate transformed divF
% for i in range(nupts):
	% for j, ex in enumerate(srcex):
	    tdivtconf[${i}][${j}] = -rcpdjac*tdivtconf[${i}][${j}] + ${ex};
	% endfor
% endfor

if (shockcell == 1){

	// Calculate solpred = forward Euler approximation of u[shockvar] 
	% for i in range(nupts):
		% for k in range(nvars):
			solpred[${i}][${k}] = u[${i}][${k}] + ${dt}*tdivtconf[${i}][${k}];
		% endfor
	% endfor

	// Calculate solpredmodes = modal coefficients of solpred
	% for i in range(nupts):
		% for k in range(nvars):
			tmp = ${' + '.join('{jx}*solpred[{j}][{k}]'.format(j=j, jx=jx, k=k) for j, jx in enumerate(invvdm[i]) if jx != 0)};
			solpredmodes[${i}][${k}] = tmp;
		% endfor
	% endfor


	// Calculate ratio of solution magnitudes of shockvar to target ratio
	// ubdegs = total order term (index offset + 1 since ignoring mean mode)
	// Target ratio = k^-4, value starts at 2 since mode ratios start at 1
	% for i in range(nupts-1):
		targetratio = pow(${ubdegs[i+1] + 1.0}, ${kexp});
		magratios[${i}] = pow(solpredmodes[${i+1}][${svar}]/solpredmodes[0][${svar}], 2.0)/targetratio;
	% endfor

	% if filtermethod == 'pointwise':
		// Compute filtercoeffs = pointwise filter coefficients
		filtercoeffs[0] = 1.0; // Leave mean mode unfiltered
		% for i in range(nupts-1):
				filtercoeffs[${i+1}] = fmax(0.0, 1.0/magratios[${i}]);
			}
			else {
				filtercoeffs[${i+1}] = 1.0;
			}
		% endfor

	% elif filtermethod == 'linear':
		// Compute filtercoeffs = pointwise filter coefficients
		fpdtype_t maxslope = 0.0, slope;
		% for i in range(nupts-1):			
			if (magratios[${i}] > 1.0){
				slope = fmax(0.0, 1.0 - 1.0/magratios[${i}])/${ubdegs[i+1]};
				if (slope > maxslope){ 
					maxslope = slope;
				}
			}
		% endfor

		filtercoeffs[0] = 1.0; // Leave mean mode unfiltered
		% for i in range(nupts-1):
			filtercoeffs[${i+1}] = fmax(0.0, 1.0 - maxslope*${ubdegs[i+1]});
		% endfor

	% endif

	//printf("%f\n", magratios[2]);
	// Compute filteredsolmodes = filtered solution modes
	% for i in range(nupts):
		% for k in range(nvars):
			filteredsolmodes[${i}][${k}] = filtercoeffs[${i}]*solpredmodes[${i}][${k}];
		% endfor
	% endfor
	//printf("%f\n", filteredsolmodes[0][0] );

	// Compute filteredsol = nodal version of filtered solution modes and then compute filtered divF using forward Euler approximation
	% for i in range(nupts):
		% for k in range(nvars):
			tmp = ${' + '.join('{jx}*filteredsolmodes[{j}][{k}]'.format(j=j, jx=jx, k=k) for j, jx in enumerate(vdm[i]) if jx != 0)};
			tdivtconf[${i}][${k}] = (tmp - u[${i}][${k}])/${dt};
		% endfor
	% endfor


}




</%pyfr:kernel>
