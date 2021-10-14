# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='zero_interior_pressure' ndim='2'
              u='inout fpdtype_t[${str(nvars)}]'>

    
    u[${ndims}] = 0.0; 

</%pyfr:kernel>
