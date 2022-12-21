# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% inf = 1e20 %>
<%pyfr:macro name='compute_entropy' params='u, d, p, e'>
    d = u[0];
    fpdtype_t rcpd = 1.0/d;
    fpdtype_t E = u[${nvars - 1}];

    fpdtype_t rhov2 = ${pyfr.dot('u[{i}+1]', i=ndims)};

    % if ndims == 2:
    fpdtype_t BdotB2 = 0.5*(u[3]*u[3] + u[4]*u[4]);
    % elif ndims == 3:
    fpdtype_t BdotB2 = 0.5*(u[4]*u[4] + u[5]*u[5] + u[6]*u[6]);
    % endif

    // Compute the pressure
    % if eos == 'ideal':
    p = ${c['gamma'] - 1}*(E - 0.5*rcpd*(rhov2) - BdotB2);
    // Compute specific physical entropy
    e = ((d > 0) && (p > 0)) ? p*pow(rcpd, ${c['gamma']}) : ${inf};
    % elif eos == 'adiabatic':
    p = ${c['kappa']}*pow(d, ${c['gamma']});
    // Compute specific physical entropy
    e = ((d > 0) && (p > 0)) ? d*log(p*pow(rcpd, ${c['gamma']})) : ${inf};
    % endif

</%pyfr:macro>
