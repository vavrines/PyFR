# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.smats'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>
<%include file='pyfr.solvers.mpnavstokes.kernels.flux'/>

<%pyfr:kernel name='tfluxlin' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              artvisc='in broadcast-col fpdtype_t'
              grad='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
              f='out fpdtype_t[${str(ndims)}][${str(nvars)}]'
              verts='in broadcast-col fpdtype_t[${str(nverts)}][${str(ndims)}]'
              upts='in broadcast-row fpdtype_t[${str(ndims)}]'>
    // Compute the flux (F = Fi + Fv)
    fpdtype_t ftemp[${ndims}][${nvars}];
    fpdtype_t a[${nspec}], d, p, v[${ndims}];
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp',  'a', 'd', 'p', 'v')};
    ${pyfr.expand('viscous_flux_add', 'u', 'grad', 'ftemp')};

    // Compute the S matrices
    fpdtype_t smats[${ndims}][${ndims}], djac;
    ${pyfr.expand('calc_smats_detj', 'verts', 'upts', 'smats', 'djac')};

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ${' + '.join(f'smats[{i}][{k}]*ftemp[{k}][{j}]'
                                 for k in range(ndims))};
% endfor
</%pyfr:kernel>
