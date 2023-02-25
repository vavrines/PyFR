# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='rsolve' params='fl, fr, n, nF, u'>
    fpdtype_t un;
    for (int i = 0; i < ${nvars}; i++) {        
        un = ${pyfr.dot('u[i][{j}]', 'n[{j}]', j=ndims)};
        nF[i] = un > 0.0 ? un*fl[i] : un*fr[i];
    }
</%pyfr:macro>
