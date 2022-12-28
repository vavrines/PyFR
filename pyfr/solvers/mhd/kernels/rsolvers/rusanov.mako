<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mhd.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;

    fpdtype_t rl = ul[0];
    fpdtype_t rr = ur[0];
    ${pyfr.expand('ideal_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('ideal_flux', 'ur', 'fr', 'pr', 'vr')};

    // Estimate the sound speed squared
    fpdtype_t c2 = ${c['gamma']}*(pl + pr)/(ul[0] + ur[0]);

    // Get the normal velocities
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    % if ndims == 2:
    fpdtype_t Bl[${ndims}] = {ul[3], ul[4]};
    fpdtype_t Br[${ndims}] = {ur[3], ur[4]};
    % elif ndims == 3:
    fpdtype_t Bl[${ndims}] = {ul[4], ul[5], ul[6]};
    fpdtype_t Br[${ndims}] = {ur[4], ur[5], ur[6]};
    % endif

    fpdtype_t BdotBl = ${pyfr.dot('Bl[{i}]', 'Bl[{i}]', i=ndims)};
    fpdtype_t BdotBr = ${pyfr.dot('Br[{i}]', 'Br[{i}]', i=ndims)};

    fpdtype_t Bldotn = ${pyfr.dot('Bl[{i}]', 'n[{i}]', i=ndims)};
    fpdtype_t Brdotn = ${pyfr.dot('Br[{i}]', 'n[{i}]', i=ndims)};

    // Compute common factors
    fpdtype_t gl2 = (${c['gamma']}*pl + BdotBl)/rl;
    fpdtype_t gr2 = (${c['gamma']}*pr + BdotBr)/rr;

    // Compute max magneto-acoustic wavespeeds
    fpdtype_t cl2 = 0.5*(gl2 + sqrt(gl2*gl2 - ${4*c['gamma']}*pl*Bldotn*Bldotn/(rl*rl)));
    fpdtype_t cr2 = 0.5*(gr2 + sqrt(gr2*gr2 - ${4*c['gamma']}*pr*Brdotn*Brdotn/(rr*rr)));
    fpdtype_t c = max(sqrt(cl2), sqrt(cr2));

    // Get max wavespeed
    fpdtype_t lambda = max(fabs(nvl) + sqrt(cl2), fabs(nvr) + sqrt(cr2));
    
    // Output
% for i in range(nvars):
    % if i == 2*ndims + 1 and divmethod == 'global':
    // Centered approximation for computing divB
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))}); 
    % elif i == 2*ndims + 1 and divmethod == 'local':
    // Interior approximation for computing divB (no correction term)
    nf[${i}] = ${' + '.join('n[{j}]*(fl[{j}][{i}])'.format(i=i, j=j) for j in range(ndims))}; 
    % else:
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + 0.5*lambda*(ul[${i}] - ur[${i}]);
    % endif
% endfor
</%pyfr:macro>