# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='correct_pressure' ndim='1'
              uoutb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uinb='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ucpy='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ufpts='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'>


    // SET UPWIND COMMON SOLUTION TO CENTERED INSTEAD

    // Calculate mesh spacing (assuming constant J: elements are cartesian)
    fpdtype_t h = pow(1.0/rcpdjac[0], ${1.0/ndims});

    // Set divergence
    fpdtype_t divu[${nupts}];
    % for i in range(nupts):
        divu[${i}] = uoutb[${i}][${ndims}];
    % endfor

    // Calculate Laplacian of pressure from interface contributions
    fpdtype_t intp[${nupts}];
    % for i in range(nupts):
        intp[${i}] = (${' + '.join('{jx}*ufpts[{j}][{var}]'.format(j=j, jx=jx, var=ndims) for j, jx in enumerate(FLM[i]) if jx != 0)})*${dt}/(h*h); 
    % endfor

    // Compute pressure
    fpdtype_t P[${nupts}];

    // Calculate pressure Poisson equation
    % for i in range(nupts):
        P[${i}] = (${' + '.join('{jx}*(divu[{j}] - intp[{j}])'.format(j=j, jx=jx, dt=dt) for j, jx in enumerate(ILM[i]) if jx != 0)})*(h*h)/${dt}; 
    % endfor

    // Calculate uncorrected pressure gradients and add pressure gradient contributions to velocity
    % if ndims == 2:
        fpdtype_t Px[${nupts}], Py[${nupts}];
        % for i in range(nupts):
            Px[${i}] = (${' + '.join('{jx}*P[{j}]'.format(j=j, jx=jx) for j, jx in enumerate(Mx[i]) if jx != 0)})/(h);
            Py[${i}] = (${' + '.join('{jx}*P[{j}]'.format(j=j, jx=jx) for j, jx in enumerate(My[i]) if jx != 0)})/(h);

            uoutb[${i}][0] += -${dt}*(Px[${i}]);
            uoutb[${i}][1] += -${dt}*(Py[${i}]);
        % endfor
    % elif ndims == 3:
        fpdtype_t Px[${nupts}], Py[${nupts}], Pz[${nupts}];
        % for i in range(nupts):
            Px[${i}] = (${' + '.join('{jx}*P[{j}]'.format(j=j, jx=jx) for j, jx in enumerate(Mx[i]) if jx != 0)})/(h);
            Py[${i}] = (${' + '.join('{jx}*P[{j}]'.format(j=j, jx=jx) for j, jx in enumerate(My[i]) if jx != 0)})/(h);
            Pz[${i}] = (${' + '.join('{jx}*P[{j}]'.format(j=j, jx=jx) for j, jx in enumerate(Mz[i]) if jx != 0)})/(h);

            uoutb[${i}][0] += -${dt}*Px[${i}];
            uoutb[${i}][1] += -${dt}*Py[${i}];
            uoutb[${i}][2] += -${dt}*Pz[${i}];
        % endfor
    % endif

    // DEBUGGING
    //fpdtype_t tmp;
    % for i in range(nupts):
        //printf("%f\n", intp[${i}]); 
        //printf("%f\n", divu[${i}]); 
        //printf("%f\n", divu[${i}]-intp[${i}]); 
        //printf("%f\n", P[${i}]); 
        //printf("%f\n", Px[${i}]); 
        //printf("%f\n", Py[${i}]); 
        //tmp =  (${' + '.join('{jx}*uoutb[{j}][0]'.format(j=j, jx=jx) for j, jx in enumerate(Lx[i]) if jx != 0)});
        //tmp += (${' + '.join('{jx}*uoutb[{j}][1]'.format(j=j, jx=jx) for j, jx in enumerate(Ly[i]) if jx != 0)});
        //printf("%f\n", tmp); 
    % endfor
    % for i in range(nfpts):
        //printf("%f\n", ufpts[${i}][2]); 
    % endfor
    //printf("%f\n", h); 
    % for i in range(nupts):
        //printf("%f\n", intp[${i}]); 
        //printf("%f\n", divu[${i}]); 
        //printf("%f\n", divu[${i}]-intp[${i}]); 
        //printf("%f\n", P[${i}]); 
        //printf("%f\n", Px[${i}]); 
        //printf("%f\n", Py[${i}]); 
    % endfor


    // Compute divF from forward Euler approximation of next state (store in uoutb)
    % for i,j in pyfr.ndrange(nupts, nvars-1):
        uoutb[${i}][${j}] = (uoutb[${i}][${j}] - uinb[${i}][${j}])/${dt};
    % endfor

    // Set pressure for uinb, set divergence to 0 for uoutb
    % for i in range(nupts):
        uinb[${i}][${nvars-1}] = P[${i}];
        uoutb[${i}][${nvars-1}] = 0.0;
    % endfor



</%pyfr:kernel>
