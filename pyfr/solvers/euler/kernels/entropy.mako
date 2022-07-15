# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% inf = 1e20 %>
<%pyfr:macro name='compute_entropy' params='u, d, p, e'>
    d = u[0];
    fpdtype_t E = u[${nvars - 1}];

    // Compute the pressure
    p = ${c['gamma'] - 1}*(E - 0.5*(${pyfr.dot('u[{i}+1]', i=ndims)})/d);

    // Compute specific physical entropy
    e = ((d > 0) && (p > 0)) ? d*(log(p) - ${c['gamma']}*log(d)) : ${inf};
</%pyfr:macro>
