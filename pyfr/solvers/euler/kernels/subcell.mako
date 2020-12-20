# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='subcell' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              plocu='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              divf='out fpdtype_t[${str(nupts)}][${str(nvars)}]'
              scverts='in fpdtype_t[${str((order+2)**2)}][${str(ndims)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              >

<%pyfr:macro name='get_tangent' params='xl,xr,n,nmag'>
    nmag = 0;
    % for dim in range(ndims):
        n[${dim}] = xl[${dim}] - xr[${dim}];
    % endfor        

    // Normalize norm
    % for dim in range(ndims):
        nmag += pow(n[${dim}], 2);
    % endfor
    nmag = pow(nmag, 0.5);
    % for dim in range(ndims):
        n[${dim}] /= nmag; 
    % endfor
</%pyfr:macro>

<%pyfr:macro name='get_normal' params='xl,xr,t,s'>
    // Counter-clockwise rotation
    fpdtype_t tmp;
    ${pyfr.expand('get_tangent','xl', 'xr', 't', 's')}

    tmp = t[1];
    t[1] = -t[0];
    t[0] = tmp;
</%pyfr:macro>

<%pyfr:macro name='get_mass' params='x1,x2,x3,x4,area'>
    area =  0.5*(x1[0]*x2[1] + x2[0]*x3[1] + x3[0]*x4[1] + x4[0]*x1[1])
          - 0.5*(x1[1]*x2[0] + x2[1]*x3[0] + x3[1]*x4[0] + x4[1]*x1[0]);

</%pyfr:macro>


% for i,j in pyfr.ndrange(nupts, nvars):
    divf[${i}][${j}] = 0.0;
% endfor

fpdtype_t fntemp[${nvars}];
fpdtype_t uc[${nvars}], u2[${nvars}];
fpdtype_t xbl[${ndims}], xbr[${ndims}], xtl[${ndims}], xtr[${ndims}];
fpdtype_t nb[${ndims}], nr[${ndims}], nt[${ndims}], nl[${ndims}];
fpdtype_t sl, sr, st, sb, mass;


% for i in range(order+1):
    % for j in range(order+1):

        % for var in range(nvars):
            uc[${var}] = u[${i+j*(order+1)}][${var}];
        % endfor

        % for dim in range(ndims):
            xbl[${dim}] = scverts[${i+(j)*(order+1)}][${dim}];
            xbr[${dim}] = scverts[${i+1+(j)*(order+1)}][${dim}];
            xtl[${dim}] = scverts[${i+(j+1)*(order+1)}][${dim}];
            xtr[${dim}] = scverts[${i+1+(j+1)*(order+1)}][${dim}];
        % endfor

        ${pyfr.expand('get_normal','xbl', 'xbr', 'nb', 'sb')}
        ${pyfr.expand('get_normal','xbr', 'xtr', 'nr', 'sr')}
        ${pyfr.expand('get_normal','xtr', 'xtl', 'nt', 'st')}
        ${pyfr.expand('get_normal','xtl', 'xbl', 'nl', 'sl')}
        ${pyfr.expand('get_mass','xbl', 'xbr', 'xtr', 'xtl', 'mass')}



        % if i != 0:
            % for var in range(nvars):
                u2[${var}] = u[${i-1+j*(order+1)}][${var}];
            % endfor
        % else:
            % for var in range(nvars):
                u2[${var}] = uc[${var}];
            % endfor
        % endif

        ${pyfr.expand('rsolve','uc','u2','nl','fntemp')};

        % for var in range(nvars):
            divf[${i+j*(order+1)}][${var}] += fntemp[${var}]*sl/mass;
        % endfor

        % if i != order:
            % for var in range(nvars):
                u2[${var}] = u[${i+1+j*(order+1)}][${var}];
            % endfor
        % else:
            % for var in range(nvars):
                u2[${var}] = uc[${var}];
            % endfor
        % endif


        ${pyfr.expand('rsolve','uc','u2','nr','fntemp')};

        % for var in range(nvars):
            divf[${i+j*(order+1)}][${var}] += fntemp[${var}]*sr/mass;
        % endfor


        % if j != 0:
            % for var in range(nvars):
                u2[${var}] = u[${i+(j-1)*(order+1)}][${var}];
            % endfor
        % else:
            % for var in range(nvars):
                u2[${var}] = uc[${var}];
            % endfor
        % endif


        ${pyfr.expand('rsolve','uc','u2','nb','fntemp')};

        % for var in range(nvars):
            divf[${i+j*(order+1)}][${var}] += fntemp[${var}]*sb/mass;
        % endfor

        % if j != order:
            % for var in range(nvars):
                u2[${var}] = u[${i+(j+1)*(order+1)}][${var}];
            % endfor
        % else:
            % for var in range(nvars):
                u2[${var}] = uc[${var}];
            % endfor
        % endif


        ${pyfr.expand('rsolve','uc','u2','nt','fntemp')};

        % for var in range(nvars):
            divf[${i+j*(order+1)}][${var}] += fntemp[${var}]*st/mass;
        % endfor

    % endfor
% endfor

% for i,j in pyfr.ndrange(nupts, nvars):
    divf[${i}][${j}] /= -rcpdjac[${i}];
% endfor



</%pyfr:kernel>
