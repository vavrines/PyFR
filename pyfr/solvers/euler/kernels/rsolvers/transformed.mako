# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% t_tol = 0.99 %>
// Transforms to m=[1,0,0]^T, where u[0],u[end] are not modified
// See Moler and Hughes 1999
<%pyfr:macro name='transform_to' params='n,u,t'>

% if ndims == 2:

    t[0] = u[0];
    t[1] =  n[0]*u[1] + n[1]*u[2];
    t[2] = -n[1]*u[1] + n[0]*u[2];
    t[3] = u[3];
    
% elif ndims == 3:

    t[0] = u[0];
    t[4] = u[4];

    if (fabs(n[0]) < ${t_tol}){
        fpdtype_t h = 1./(1. + n[0]);

        t[1] =  n[0]*u[1] + n[1]*u[2] + n[2]*u[3];
    	t[2] = -n[1]*u[1] + (n[0] + h*n[2]*n[2])*u[2] - h*n[1]*n[2]*u[3];
    	t[3] = -n[2]*u[1] - h*n[1]*n[2]*u[2] + (n[0] + h*n[1]*n[1])*u[3];
    }
    else if (fabs(n[1]) < fabs(n[2])){
        fpdtype_t h = 1./(1. - n[1]);
	
        t[1] = n[0]*u[1] + n[1]*u[2] + n[2]*u[3];
	t[2] =  (1. - h*n[0]*n[0])*u[1] + n[0]*u[2] - h*n[0]*n[2]*u[3];
	t[3] = -h*n[0]*n[2]*u[1] + n[2]*u[2] + (1. - h*n[2]*n[2])*u[3];
    }
    else{
       fpdtype_t h = 1./(1. - n[2]);
       
       t[1] = n[0]*u[1] + n[1]*u[2] + n[2]*u[3];
       t[2] = -h*n[0]*n[1]*u[1] + (1. - h*n[1]*n[1])*u[2] + n[1]*u[3];
       t[3] =  (1. - h*n[0]*n[0])*u[1] - h*n[0]*n[1]*u[2] + n[0]*u[3];
    }

% endif
</%pyfr:macro>

// Transforms from m=[1,0,0]^T, where u[0],u[end] are not modified
// See Moler and Hughes 1999
<%pyfr:macro name='transform_from' params='n,t,u'>

% if ndims == 2:

    u[0] = t[0];
    u[1] = n[0]*t[1] - n[1]*t[2];
    u[2] = n[1]*t[1] + n[0]*t[2];
    u[3] = t[3];
    
% elif ndims == 3:

    u[0] =  t[0];
    u[4] =  t[4];

    if (fabs(n[0]) < ${t_tol}){
        fpdtype_t h = 1./(1. + n[0]);

        u[1] =  n[0]*t[1] - n[1]*t[2] - n[2]*t[3];
    	u[2] =  n[1]*t[1] + (n[0] + h*n[2]*n[2])*t[2] - h*n[1]*n[2]*t[3];
    	u[3] =  n[2]*t[1] - h*n[1]*n[2]*t[2] + (n[0] + h*n[1]*n[1])*t[3];
    }
    else if (fabs(n[1]) < fabs(n[2])){
        fpdtype_t h = 1./(1. - n[1]);
	
        u[1] = n[0]*t[1] +  (1. - h*n[0]*n[0])*t[2] - h*n[0]*n[2]*t[3];
	u[2] = n[1]*t[1] + n[0]*t[2] + n[2]*t[3];
	u[3] = n[2]*t[1] - h*n[0]*n[2]*t[2] + (1. - h*n[2]*n[2])*t[3];
    }
    else{
       fpdtype_t h = 1./(1. - n[2]);
       
       u[1] = n[0]*t[1] - h*n[0]*n[1]*t[2] + (1. - h*n[0]*n[0])*t[3];
       u[2] = n[1]*t[1] + (1. - h*n[1]*n[1])*t[2] - h*n[0]*n[1]*t[3];
       u[3] = n[2]*t[1] + n[1]*t[2] + n[0]*t[3];
    }


% endif
</%pyfr:macro>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];

    ${pyfr.expand('transform_to','n', 'ul', 'utl')};
    ${pyfr.expand('transform_to','n', 'ur', 'utr')};

    ${pyfr.expand('rsolve_t1d','utl','utr','ntf')};

    ${pyfr.expand('transform_from','n','ntf','nf')};

</%pyfr:macro>