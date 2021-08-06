# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
</%pyfr:macro>

% if ndims == 2:
    <%pyfr:macro name='inviscid_flux' params='uin, fout'>
        fout[0][0] = sin(uin[0]);
        fout[1][0] = cos(uin[0]);
    </%pyfr:macro>
% elif ndims == 3:
    <%pyfr:macro name='inviscid_flux' params='uin, fout'>
        fout[0][0] = sin(uin[0]);
        fout[1][0] = cos(uin[0]);
        fout[2][0] = 0.0;
    </%pyfr:macro>
% endif

