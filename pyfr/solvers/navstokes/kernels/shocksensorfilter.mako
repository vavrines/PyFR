# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<%pyfr:kernel name='shocksensorfilter' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              shockcell='out fpdtype_t'>

shockcell = 1;



</%pyfr:kernel>
