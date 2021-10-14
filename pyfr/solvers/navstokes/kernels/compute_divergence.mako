# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='compute_divergence' ndim='2'
              u='inout fpdtype_t[${str(nvars)}]'
              gradu='in fpdtype_t[${str(ndims)}][${str(nvars)}]'>

    % if ndims == 2:
        u[${ndims}] = gradu[0][0] + gradu[1][1];
    % elif ndims == 3:
        u[${ndims}] = gradu[0][0] + gradu[1][1] + gradu[2][2];
    % endif


</%pyfr:kernel>
