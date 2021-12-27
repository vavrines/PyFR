# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='calcminentropy' ndim='1'
              entmin_int='in fpdtype_t[${str(nfpts)}]'
              entmin='out fpdtype_t'>

fpdtype_t delta_ent = 0;
% for i in range(nfpts):
    entmin = fmin(entmin, entmin_int[${i}]);
    delta_ent += entmin_int[${i}];
% endfor

delta_ent = ${alpha}*(delta_ent/${nfpts} - entmin);
entmin -= fmax(0, delta_ent);

</%pyfr:kernel>
