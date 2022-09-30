# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

% if ndims == 2:

<% i_u = nspec + 0 %>
<% i_v = nspec + 1 %>
<% i_E = nspec + 2 %>

<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
    fpdtype_t rho = ${' + '.join('uin[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rhou = uin[${i_u}], rhov = uin[${i_v}], E = uin[${i_E}];

    fpdtype_t a[${nspec}];
% for i in range(nspec - 1):
    a[${i}] = uin[${nspec + ndims + 1 + i}];
% endfor
    a[${nspec - 1}] = 1 - (${' + '.join('a[{i}]'.format(i=i) for i in range(nspec - 1))});

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov;

    fpdtype_t rho_x = ${' + '.join('grad_uin[0][{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rho_y = ${' + '.join('grad_uin[1][{i}]'.format(i=i) for i in range(nspec))};

    // Velocity derivatives (rho*grad[u,v])
    fpdtype_t u_x = grad_uin[0][${i_u}] - u*rho_x;
    fpdtype_t u_y = grad_uin[1][${i_u}] - u*rho_y;
    fpdtype_t v_x = grad_uin[0][${i_v}] - v*rho_x;
    fpdtype_t v_y = grad_uin[1][${i_v}] - v*rho_y;

    fpdtype_t E_x = grad_uin[0][${i_E}];
    fpdtype_t E_y = grad_uin[1][${i_E}];

    fpdtype_t mu_c = ${' + '.join('a[{i}]*{mu}'.format(i=i, mu=c[f'mu{i}']) for i in range(nspec))};

    // Compute temperature derivatives (c_v*dT/d[x,y])
    fpdtype_t T_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x));
    fpdtype_t T_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y));

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*mu_c*rcprho*(u_x - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_yy = -2*mu_c*rcprho*(v_y - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_xy = -mu_c*rcprho*(v_x + u_y);

    fout[0][${i_u}] += t_xx; fout[1][${i_u}] += t_xy;
    fout[0][${i_v}] += t_xy; fout[1][${i_v}] += t_yy;

    fpdtype_t kap = ${' + '.join('a[{i}]*{kap}'.format(i=i, kap=(c[f'mu{i}']*c[f'gamma{i}']/c[f'Pr{i}'])) for i in range(nspec))};

    fout[0][${i_E}] += u*t_xx + v*t_xy - kap*T_x;
    fout[1][${i_E}] += u*t_xy + v*t_yy - kap*T_y;
</%pyfr:macro>
% elif ndims == 3:

<% i_u = nspec + 0 %>
<% i_v = nspec + 1 %>
<% i_w = nspec + 2 %>
<% i_E = nspec + 3 %>

<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
    fpdtype_t rho = ${' + '.join('uin[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rhou = uin[${i_u}], rhov = uin[${i_v}], rhow = uin[${i_w}]
    fpdtype_t E = uin[${i_E}];

    fpdtype_t a[${nspec}];
% for i in range(nspec - 1):
    a[${i}] = uin[${nspec + ndims + 1 + i}];
% endfor
    a[${nspec - 1}] = 1 - (${' + '.join('a[{i}]'.format(i=i) for i in range(nspec - 1))});

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov, w = rcprho*rhow;

    fpdtype_t rho_x = ${' + '.join('grad_uin[0][{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rho_y = ${' + '.join('grad_uin[1][{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rho_z = ${' + '.join('grad_uin[2][{i}]'.format(i=i) for i in range(nspec))};

    // Velocity derivatives (rho*grad[u,v,w])
    fpdtype_t u_x = grad_uin[0][${i_u}] - u*rho_x;
    fpdtype_t u_y = grad_uin[1][${i_u}] - u*rho_y;
    fpdtype_t u_z = grad_uin[2][${i_u}] - u*rho_z;
    fpdtype_t v_x = grad_uin[0][${i_v}] - v*rho_x;
    fpdtype_t v_y = grad_uin[1][${i_v}] - v*rho_y;
    fpdtype_t v_z = grad_uin[2][${i_v}] - v*rho_z;
    fpdtype_t w_x = grad_uin[0][${i_w}] - w*rho_x;
    fpdtype_t w_y = grad_uin[1][${i_w}] - w*rho_y;
    fpdtype_t w_z = grad_uin[2][${i_w}] - w*rho_z;

    fpdtype_t E_x = grad_uin[0][${i_E}];
    fpdtype_t E_y = grad_uin[1][${i_E}];
    fpdtype_t E_z = grad_uin[2][${i_E}];

    fpdtype_t mu_c = ${' + '.join('a[{i}]*{mu}'.format(i=i, mu=c[f'mu{i}']) for i in range(nspec))};

    // Compute temperature derivatives (c_v*dT/d[x,y,z])
    fpdtype_t T_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x + w*w_x));
    fpdtype_t T_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y + w*w_y));
    fpdtype_t T_z = rcprho*(E_z - (rcprho*rho_z*E + u*u_z + v*v_z + w*w_z));

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*mu_c*rcprho*(u_x - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_yy = -2*mu_c*rcprho*(v_y - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_zz = -2*mu_c*rcprho*(w_z - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_xy = -mu_c*rcprho*(v_x + u_y);
    fpdtype_t t_xz = -mu_c*rcprho*(u_z + w_x);
    fpdtype_t t_yz = -mu_c*rcprho*(w_y + v_z);

    fout[0][${i_u}] += t_xx; fout[1][${i_u}] += t_xy; fout[2][${i_u}] += t_xz;
    fout[0][${i_v}] += t_xy; fout[1][${i_v}] += t_yy; fout[2][${i_v}] += t_yz;
    fout[0][${i_w}] += t_xz; fout[1][${i_w}] += t_yz; fout[2][${i_w}] += t_zz;

    fpdtype_t kap = ${' + '.join('a[{i}]*{kap}'.format(i=i, kap=(c[f'mu{i}']*c[f'gamma{i}']/c[f'Pr{i}'])) for i in range(nspec))};

    fout[0][${i_E}] += u*t_xx + v*t_xy + w*t_xz - kap*T_x;
    fout[1][${i_E}] += u*t_xy + v*t_yy + w*t_yz - kap*T_y;
    fout[2][${i_E}] += u*t_xz + v*t_yz + w*t_zz - kap*T_z;
</%pyfr:macro>
% endif
