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
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr, pgl, pgr;

    fpdtype_t rl = ul[0];
    fpdtype_t rr = ur[0];
    ${pyfr.expand('ideal_flux', 'ul', 'fl', 'pgl', 'vl')};
    ${pyfr.expand('ideal_flux', 'ur', 'fr', 'pgr', 'vr')};

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
    fpdtype_t gl2 = (${c['gamma']}*pgl + BdotBl)/rl;
    fpdtype_t gr2 = (${c['gamma']}*pgr + BdotBr)/rr;

    // Compute max magneto-acoustic wavespeeds
    fpdtype_t cl2 = 0.5*(gl2 + sqrt(gl2*gl2 - ${4*c['gamma']}*pgl*Bl[0]*Bl[0]/(rl*rl)));
    fpdtype_t cr2 = 0.5*(gr2 + sqrt(gr2*gr2 - ${4*c['gamma']}*pgr*Br[0]*Br[0]/(rr*rr)));
    fpdtype_t c = max(sqrt(cl2), sqrt(cr2));

    // Get the normal left and right velocities
    fpdtype_t nvl = vl[0];
    fpdtype_t nvr = vr[0];

    // Estimate the left and right wave speed, sl and sr
    fpdtype_t sl = nvl - c;
    fpdtype_t sr = nvr + c;
    
    if (sl >= 0) {
        % for i in range(nvars):
        nf[${i}] = fl[0][${i}];
        % endfor
    }
    else if (sl <= 0 && sr >= 0) {
        % for i in range(nvars):
        nf[${i}] = (sr*fl[0][${i}] - sl*fr[0][${i}] + sl*sr*(ur[${i}] - ul[${i}]))/(sr - sl);
        % endfor
    }
    else {
        % for i in range(nvars):
        nf[${i}] = fr[0][${i}];
        % endfor
    }

    % if divmethod == 'global':
    // Centered approximation for computing divB
    nf[${2*ndims + 1}] = 0.5*(fl[0][${2*ndims + 1}] + fr[0][${2*ndims + 1}]);
    % elif divmethod == 'local':
    // Interior approximation for computing divB (no correction term)
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

</%pyfr:macro>