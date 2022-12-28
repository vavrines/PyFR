<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mhd.kernels.flux'/>

<% t_tol = 0.99 %>

// Transforms to m=[1,0,0]^T
// See Moler and Hughes 1999
<%pyfr:macro name='transform_to' params='n, u, t, o1, o2'>
% if ndims == 2:
    t[o1 + 0] =  n[0]*u[o1 + 0] + n[1]*u[o1 + 1];
    t[o1 + 1] = -n[1]*u[o1 + 0] + n[0]*u[o1 + 1];
    t[o2 + 0] =  n[0]*u[o2 + 0] + n[1]*u[o2 + 1];
    t[o2 + 1] = -n[1]*u[o2 + 0] + n[0]*u[o2 + 1];
% elif ndims == 3:
    if (fabs(n[0]) < ${t_tol})
    {
        fpdtype_t h = 1/(1 + n[0]);
        t[o1 + 0] =  n[0]*u[o1 + 0] + n[1]*u[o1 + 1] + n[2]*u[o1 + 2];
        t[o1 + 1] = -n[1]*u[o1 + 0] + (n[0] + h*n[2]*n[2])*u[o1 + 1] - h*n[1]*n[2]*u[o1 + 2];
        t[o1 + 2] = -n[2]*u[o1 + 0] - h*n[1]*n[2]*u[o1 + 1] + (n[0] + h*n[1]*n[1])*u[o1 + 2];
        t[o2 + 0] =  n[0]*u[o2 + 0] + n[1]*u[o2 + 1] + n[2]*u[o2 + 2];
        t[o2 + 1] = -n[1]*u[o2 + 0] + (n[0] + h*n[2]*n[2])*u[o2 + 1] - h*n[1]*n[2]*u[o2 + 2];
        t[o2 + 2] = -n[2]*u[o2 + 0] - h*n[1]*n[2]*u[o2 + 1] + (n[0] + h*n[1]*n[1])*u[o2 + 2];
    }
    else if (fabs(n[1]) < fabs(n[2]))
    {
        fpdtype_t h = 1/(1 - n[1]);
        t[o1 + 0] = n[0]*u[o1 + 0] + n[1]*u[o1 + 1] + n[2]*u[o1 + 2];
        t[o1 + 1] =  (1 - h*n[0]*n[0])*u[o1 + 0] + n[0]*u[o1 + 1] - h*n[0]*n[2]*u[o1 + 2];
        t[o1 + 2] = -h*n[0]*n[2]*u[o1 + 0] + n[2]*u[o1 + 1] + (1 - h*n[2]*n[2])*u[o1 + 2];
        t[o2 + 0] = n[0]*u[o2 + 0] + n[1]*u[o2 + 1] + n[2]*u[o2 + 2];
        t[o2 + 1] =  (1 - h*n[0]*n[0])*u[o2 + 0] + n[0]*u[o2 + 1] - h*n[0]*n[2]*u[o2 + 2];
        t[o2 + 2] = -h*n[0]*n[2]*u[o2 + 0] + n[2]*u[o2 + 1] + (1 - h*n[2]*n[2])*u[o2 + 2];
    }
    else
    {
        fpdtype_t h = 1/(1 - n[2]);
        t[o1 + 0] = n[0]*u[o1 + 0] + n[1]*u[o1 + 1] + n[2]*u[o1 + 2];
        t[o1 + 1] = -h*n[0]*n[1]*u[o1 + 0] + (1 - h*n[1]*n[1])*u[o1 + 1] + n[1]*u[o1 + 2];
        t[o1 + 2] =  (1 - h*n[0]*n[0])*u[o1 + 0] - h*n[0]*n[1]*u[o1 + 1] + n[0]*u[o1 + 2];
        t[o2 + 0] = n[0]*u[o2 + 0] + n[1]*u[o2 + 1] + n[2]*u[o2 + 2];
        t[o2 + 1] = -h*n[0]*n[1]*u[o2 + 0] + (1 - h*n[1]*n[1])*u[o2 + 1] + n[1]*u[o2 + 2];
        t[o2 + 2] =  (1 - h*n[0]*n[0])*u[o2 + 0] - h*n[0]*n[1]*u[o2 + 1] + n[0]*u[o2 + 2];
    }
% endif
</%pyfr:macro>

// Transforms from m=[1,0,0]^T
<%pyfr:macro name='transform_from' params='n, t, u, o1, o2'>
% if ndims == 2:
    u[o1 + 0] = n[0]*t[o1 + 0] - n[1]*t[o1 + 1];
    u[o1 + 1] = n[1]*t[o1 + 0] + n[0]*t[o1 + 1];
    u[o2 + 0] = n[0]*t[o2 + 0] - n[1]*t[o2 + 1];
    u[o2 + 1] = n[1]*t[o2 + 0] + n[0]*t[o2 + 1];
% elif ndims == 3:
    u[0] = t[0];
    if (fabs(n[0]) < ${t_tol})
    {
        fpdtype_t h = 1/(1 + n[0]);
        u[o1 + 0] =  n[0]*t[o1 + 0] - n[1]*t[o1 + 1] - n[2]*t[o1 + 2];
        u[o1 + 1] =  n[1]*t[o1 + 0] + (n[0] + h*n[2]*n[2])*t[o1 + 1] - h*n[1]*n[2]*t[o1 + 2];
        u[o1 + 2] =  n[2]*t[o1 + 0] - h*n[1]*n[2]*t[o1 + 1] + (n[0] + h*n[1]*n[1])*t[o1 + 2];
        u[o2 + 0] =  n[0]*t[o2 + 0] - n[1]*t[o2 + 1] - n[2]*t[o2 + 2];
        u[o2 + 1] =  n[1]*t[o2 + 0] + (n[0] + h*n[2]*n[2])*t[o2 + 1] - h*n[1]*n[2]*t[o2 + 2];
        u[o2 + 2] =  n[2]*t[o2 + 0] - h*n[1]*n[2]*t[o2 + 1] + (n[0] + h*n[1]*n[1])*t[o2 + 2];
    }
    else if (fabs(n[1]) < fabs(n[2]))
    {
        fpdtype_t h = 1/(1 - n[1]);
        u[o1 + 0] = n[0]*t[o1 + 0] +  (1 - h*n[0]*n[0])*t[o1 + 1] - h*n[0]*n[2]*t[o1 + 2];
        u[o1 + 1] = n[1]*t[o1 + 0] + n[0]*t[o1 + 1] + n[2]*t[o1 + 2];
        u[o1 + 2] = n[2]*t[o1 + 0] - h*n[0]*n[2]*t[o1 + 1] + (1 - h*n[2]*n[2])*t[o1 + 2];
        u[o2 + 0] = n[0]*t[o2 + 0] +  (1 - h*n[0]*n[0])*t[o2 + 1] - h*n[0]*n[2]*t[o2 + 2];
        u[o2 + 1] = n[1]*t[o2 + 0] + n[0]*t[o2 + 1] + n[2]*t[o2 + 2];
        u[o2 + 2] = n[2]*t[o2 + 0] - h*n[0]*n[2]*t[o2 + 1] + (1 - h*n[2]*n[2])*t[o2 + 2];
    }
    else
    {
        fpdtype_t h = 1/(1 - n[2]);
        u[o1 + 0] = n[0]*t[o1 + 0] - h*n[0]*n[1]*t[o1 + 1] + (1 - h*n[0]*n[0])*t[o1 + 2];
        u[o1 + 1] = n[1]*t[o1 + 0] + (1 - h*n[1]*n[1])*t[o1 + 1] - h*n[0]*n[1]*t[o1 + 2];
        u[o1 + 2] = n[2]*t[o1 + 0] + n[1]*t[o1 + 1] + n[0]*t[o1 + 2];
        u[o2 + 0] = n[0]*t[o2 + 0] - h*n[0]*n[1]*t[o2 + 1] + (1 - h*n[0]*n[0])*t[o2 + 2];
        u[o2 + 1] = n[1]*t[o2 + 0] + (1 - h*n[1]*n[1])*t[o2 + 1] - h*n[0]*n[1]*t[o2 + 2];
        u[o2 + 2] = n[2]*t[o2 + 0] + n[1]*t[o2 + 1] + n[0]*t[o2 + 2];
    }
% endif
</%pyfr:macro>

<%pyfr:macro name='rsolve_t1d' params='ul, ur, nf'>
    // Compute the left and right fluxes + velocities and ptressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t ptl, ptr, pl, pr;

    fpdtype_t rl = ul[0];
    fpdtype_t rr = ur[0];
    ${pyfr.expand('ideal_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('ideal_flux', 'ur', 'fr', 'pr', 'vr')};

    % if ndims == 2:
    fpdtype_t Bl[${ndims}] = {ul[3], ul[4]};
    fpdtype_t Br[${ndims}] = {ur[3], ur[4]};
    % elif ndims == 3:
    fpdtype_t Bl[${ndims}] = {ul[4], ul[5], ul[6]};
    fpdtype_t Br[${ndims}] = {ur[4], ur[5], ur[6]};
    % endif

    fpdtype_t BdotBl = ${pyfr.dot('Bl[{i}]', 'Bl[{i}]', i=ndims)};
    fpdtype_t BdotBr = ${pyfr.dot('Br[{i}]', 'Br[{i}]', i=ndims)};

    fpdtype_t Bdotvl = ${pyfr.dot('Bl[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t Bdotvr = ${pyfr.dot('Br[{i}]', 'vr[{i}]', i=ndims)};

    // Compute common factors
    fpdtype_t gl2 = (${c['gamma']}*pl + BdotBl)/rl;
    fpdtype_t gr2 = (${c['gamma']}*pr + BdotBr)/rr;

    // Compute max magneto-acoustic wavespeeds
    fpdtype_t cl = sqrt(0.5*(gl2 + sqrt(gl2*gl2 - ${4*c['gamma']}*pl*Bl[0]*Bl[0]/(rl*rl))));
    fpdtype_t cr = sqrt(0.5*(gr2 + sqrt(gr2*gr2 - ${4*c['gamma']}*pr*Br[0]*Br[0]/(rr*rr))));

    // Get the normal left and right velocities
    fpdtype_t nvl = vl[0];
    fpdtype_t nvr = vr[0];

    // Estimate the left and right wave speed, sl and sr
    fpdtype_t sl = min(nvl - cl, nvr - cr);
    fpdtype_t sr = max(nvl + cl, nvr + cr);

    ptl = pl + 0.5*BdotBl;
    ptr = pr + 0.5*BdotBr;
    fpdtype_t s_star = (rr*nvr*(sr - nvr) - rl*nvl*(sl - nvl) + ptl - ptr - Bl[0]*Bl[0] + Br[0]*Br[0]) / 
                       (rr*(sr - nvr) - rl*(sl - nvl));


    fpdtype_t u_hll[${nvars}];

    if (sl >= 0) {
        % for i in range(nvars):
        u_hll[${i}] = ul[${i}];
        % endfor
    }
    else if (sl < 0 && sr >= 0) {
        % for i in range(nvars):
        u_hll[${i}] = (sr*ur[${i}] - sl*ul[${i}] - fr[0][${i}] + fl[0][${i}])/(sr - sl);
        % endfor
    }
    else {
        % for i in range(nvars):
        u_hll[${i}] = ur[${i}];
        % endfor
    }

    fpdtype_t B_hll[${ndims}], v_hll[${ndims}];
    % for i in range(ndims):
    B_hll[${i}] = u_hll[${i+ndims+1}];
    v_hll[${i}] = u_hll[${i+1}]/u_hll[0];
    % endfor

    fpdtype_t Bdotv_hll = ${pyfr.dot('B_hll[{i}]', 'v_hll[{i}]', i=ndims)};
    
    if (sl >= 0) {
        % for i in range(nvars):
        nf[${i}] = fl[0][${i}];
        % endfor
    }
    else if (sl < 0 && s_star >= 0) {
        fpdtype_t ulstar[${nvars}];
        fpdtype_t p_star = rl*(sl - nvl)*(s_star - nvl) + ptl - Bl[0]*Bl[0] + B_hll[0]*B_hll[0];

        ulstar[0] = rl*(sl - nvl)/(sl - s_star);
        ulstar[1] = ulstar[0]*s_star;
        ulstar[2] = ul[2]*(sl - nvl)/(sl - s_star) - (B_hll[0]*B_hll[1] - Bl[0]*Bl[1])/(sl - s_star);

        % if ndims == 2:
        ulstar[3] = B_hll[0];
        ulstar[4] = B_hll[1];
        ulstar[5] = 0.0;
        ulstar[6] = ul[6]*(sl - nvl)/(sl - s_star) + (p_star*s_star - ptl*nvl - (B_hll[0]*Bdotv_hll - Bl[0]*Bdotvl))/(sl - s_star);
        % elif ndims == 3:
        ulstar[3] = ul[3]*(sl - nvl)/(sl - s_star) - (B_hll[0]*B_hll[2] - Bl[0]*Bl[2])/(sl - s_star);
        ulstar[4] = B_hll[0];
        ulstar[5] = B_hll[1];
        ulstar[6] = B_hll[2];
        ulstar[7] = 0.0;
        ulstar[8] = ul[8]*(sl - nvl)/(sl - s_star) + (p_star*s_star - ptl*nvl - (B_hll[0]*Bdotv_hll - Bl[0]*Bdotvl))/(sl - s_star);
        % endif
        
        % for i in range(nvars):
        nf[${i}] = fl[0][${i}] + sl*(ulstar[${i}] - ul[${i}]);
        % endfor
    }
    else if (s_star < 0 && sr >= 0) {
        fpdtype_t urstar[${nvars}];
        fpdtype_t p_star = rr*(sr - nvr)*(s_star - nvr) + ptr - Br[0]*Br[0] + B_hll[0]*B_hll[0];

        urstar[0] = rr*(sr - nvr)/(sr - s_star);
        urstar[1] = urstar[0]*s_star;
        urstar[2] = ur[2]*(sr - nvr)/(sr - s_star) - (B_hll[0]*B_hll[1] - Br[0]*Br[1])/(sr - s_star);
        
        % if ndims == 2:
        urstar[3] = B_hll[0];
        urstar[4] = B_hll[1];
        urstar[5] = 0.0;
        urstar[6] = ur[6]*(sr - nvr)/(sr - s_star) + (p_star*s_star - ptr*nvr - (B_hll[0]*Bdotv_hll - Br[0]*Bdotvr))/(sr - s_star);
        % elif ndims == 3:
        urstar[3] = ur[3]*(sr - nvr)/(sr - s_star) - (B_hll[0]*B_hll[2] - Br[0]*Br[2])/(sr - s_star);
        urstar[4] = B_hll[0];
        urstar[5] = B_hll[1];
        urstar[6] = B_hll[2];
        urstar[7] = 0.0;
        urstar[8] = ur[8]*(sr - nvr)/(sr - s_star) + (p_star*s_star - ptr*nvr - (B_hll[0]*Bdotv_hll - Br[0]*Bdotvr))/(sr - s_star);
        % endif
        
        % for i in range(nvars):
        nf[${i}] = fr[0][${i}] + sr*(urstar[${i}] - ur[${i}]);
        % endfor
    }
    else {
        % for i in range(nvars):
        nf[${i}] = fr[0][${i}];
        % endfor
    }

    % if divmethod == 'global':
    // Centered apptroximation for computing divB
    nf[${2*ndims + 1}] = 0.5*(fl[0][${2*ndims + 1}] + fr[0][${2*ndims + 1}]);
    % elif divmethod == 'local':
    // Interior apptroximation for computing divB (no correction term)
    nf[${2*ndims + 1}] = fl[0][${2*ndims + 1}];
    % endif
</%pyfr:macro>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];
    % for i in range(nvars):
    utl[${i}] = ul[${i}];
    utr[${i}] = ur[${i}];
    % endfor

    % if ndims == 2:
    ${pyfr.expand('transform_to','n', 'ul', 'utl','1','3')};
    ${pyfr.expand('transform_to','n', 'ur', 'utr','1','3')};
    % elif ndims == 3:
    ${pyfr.expand('transform_to','n', 'ul', 'utl','1','4')};
    ${pyfr.expand('transform_to','n', 'ur', 'utr','1','4')};
    % endif:

    ${pyfr.expand('rsolve_t1d','utl','utr','ntf')};
    
    % for i in range(nvars):
    nf[${i}] = ntf[${i}];
    % endfor

    % if ndims == 2:
    ${pyfr.expand('transform_from','n','ntf','nf','1','3')};
    % elif ndims == 3:
    ${pyfr.expand('transform_from','n','ntf','nf','1','4')};
    % endif

    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;

    ${pyfr.expand('ideal_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('ideal_flux', 'ur', 'fr', 'pr', 'vr')};

    
    nf[${2*ndims + 1}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=2*ndims + 1, j=j) for j in range(ndims))}); 

</%pyfr:macro>