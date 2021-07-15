# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% se0 = math.log10(c['s0']) %>

<%pyfr:kernel name='shocksensor' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              revvisc='out fpdtype_t[${str(nvars)}]'>

fpdtype_t mu = ${c['max-artvisc']};

$ for i in range(nvars):
    revvisc[${i}] = mu;
% endfor


</%pyfr:kernel>
