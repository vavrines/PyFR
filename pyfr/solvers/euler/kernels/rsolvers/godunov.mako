# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<% gmrtg = (c['gamma']-1.0)/(2.0*c['gamma']) %>
<% gprtg = (c['gamma']+1.0)/(2.0*c['gamma']) %>
<% tgrgm = (2.0*c['gamma'])/(c['gamma']-1.0) %>
<% tgrgp = (2.0*c['gamma'])/(c['gamma']+1.0) %>
<% trgm = 2.0/(c['gamma']-1.0) %>
<% trgp = 2.0/(c['gamma']+1.0) %>
<% gmrgp = (c['gamma']-1.0)/(c['gamma']+1.0) %>
<% hgm = 0.5*(c['gamma']-1.0) %>
<% rgm = 1./(c['gamma']-1.0) %>
<% gamma = c['gamma'] %>



// Initial guess for pressure
<%pyfr:macro name='init_p' params='rl,ul,pl,cl,rr,ur,pr,cr,p0'>
    fpdtype_t bpv = max(0.,0.5*(pl + pr) + 0.125*(ul - ur)*(rl + rr)*(cl + cr));
    fpdtype_t pmin = min(pl,pr);
    fpdtype_t pmax = max(pl,pr);
    fpdtype_t rpmax = pmax/pmin;

    if ((rpmax <= 2.) && (pmin <= bpv) && (bpv <= pmax)) {
       p0 = bpv;
    }
    elseif(bpv < pmin) {
    	// Two-rarefaction Riemann solve
        fpdtype_t pre = pow(pl/pr,${gmrtg});
	fpdtype_t um  = (pre*ul/cl + ur/cr + ${trgm}*(pq - 1.0))/(pre/cl + 1.0/cr);
			    
        fpdtype_t ptl = 1.0 + ${hgm}*(um - ul)/cl;
        fpdtype_t ptr = 1.0 + ${hgm}*(um - ur)/cr;
	
        p0 = 0.5*(pl*pow(ptl,${tgrgm}) + pr*pow(ptr,${tgrgm}));
    }
    else{
        // Two-shock Riemann solve	
        fpdtype_t gl = sqrt((${trgp}/rl)/(${gmrgp}*pl + bpv));
        fpdtype_t gr = sqrt((${trgp}/rr)/(${gmrgp}*pr + bpv));
	p0 = (gl*pl + gr*pr - (ur - ul))/(gt + gr);
    }
</%pyfr:macro>



// Star Flux, assuming covolume = 0. See Toro 2009 Eq.(4.86-4.87)
<%pyfr:macro name='star_flux' params='p, ps, rs, cs, f, fd'>
    if (p <= ps){
       fpdtype_t pr = p/ps;
       f  = ${trgm}*cs*(pow(prat,${gmrtg})-1.0);
       fd = pow(prat,-${gmrgp})/(rs*cs);	
    }   
    else{
       fpdtype_t as = ${trgp}/rs;
       fpdtype_t bs = ps*${gmrgp};
       fpdtype_t sapb = sqrt(as/(p+bs));
       f  = (p-ps)*sapb;
       fd = (1. - 0.5*(p-ps)/(p + bs))*sapb;
    }
</%pyfr macro>



// Primitive to inviscid flux
<%pyfr:macro name='inviscid_flux' params='w, f'>
    fpdtype_t invrho = 1.0/w[0], E;

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = w[0]*w[${i + 1}];
% endfor

    // Compute the Energy
    E = p*${rgm} + 0.5*w[0]*${pyfr.dot('w[{i+1}]', i=ndims)})

    // Density and energy fluxes
% for i in range(ndims):
    f[${i}][0] = rhov[${i}];
    f[${i}][${nvars - 1}] = (E + w[${nvars-1}])*w[${i+1}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = rhov[${i}]*w[${j+1}]${' + w[${nvars-1}]' if i == j else ''};
% endfor
</%pyfr:macro>



// Reconstruct velocity vector
<%pyfr:macro name='recvel' params='n,v0,nv0,nv1,v1'>
%for i in range(ndims):
    v1[${i}] = v0[${i}] - (nv0 - nv1)*n[${i}];
% endfor
</%pyfr:macro>



// Reconstruct primitive vector
<%pyfr:macro name='recprim' params='r,v,p,w'>
    w[0] = r;
%for i in range(ndims):
    w[${i+1}] = v[${i}];
% endfor
    w[${nvars-1}] = p;
</%pyfr:macro>


// Exact solve solution decision tree
<% switch = 0.0 %>
<%pyfr:macro name='riemann_desicion' params='n,rl,vl,pl,cl,nvl,rr,vr,pr,cr,nvr,us,p0,w0'>
    fpdtype_t vs[${ndims}];
    if (${switch} <= us){
        if (p0 <= pl){
            fpdtype_t ll = nvl - cl;
            if (${switch} <= ll){
	        ${pyfr.expand('recprim', 'rl','vl','pl','w0')};
	    }
	    else{
	        fpdtype_t cml = cl*pow(p0/pl,${gmrtg});
                fpdtype_t ls = us - cml;
		if (${switch} > ls){
		    fpdtype_t rs = rl*pow(p0/pl, 1./${gamma});
		    ${pyfr.expand('recvel','n','vl','nvl','us','vs')};
		    ${pyfr.expand('recprim','rs','vs','p0','w0')};
                }
		else{
		    fpdtype_t c = ${trgp}*(cl + ${hgm}*(nvl - ${switch}));
		    fpdtype_t rs = rl*pow(c/cl,${trgm});
		    fpdtype_t uc = ${trgp}*(cl + ${hgm}*nvl + ${switch});
		    fpdtype_t ps = pl*pow(c/cl,${tgrgm});
		    ${pyfr.expand('recvel','n','vl','nvl','uc','vs')};
		    ${pyfr.expand('recprim','rs','vs','ps','w0')};
                }
            }
        }
        else{
	    fpdtype_t p0p = p0/pl;
	    fpdtype_t sl = nvl - cl*sqrt(${gprtg}*p0p + ${gmrtg});
	    if (${switch} <= sl){
	        ${pyfr.expand('recprim','rl','vl','pl','w0')};
	    }
	    else{
	        fpdtype_t rs = rl*(p0p + ${gmrgp})/(p0p*${gmrgp} + 1.);
		${pyfr.expand('recvel','n','vl','nvl','us','vs')};
		${pyfr.expand('recprim','rs','vs','p0','w0')};
	    }
        }
    }
    else{
        if (p0 > pr){
	    fpdtype_t p0p = p0/pr;
	    fpdtype_t sr = nvr + cr*sqrt(${gprtg}*p0p + ${gmrtg});
	    if (${switch} >= sr) {
	        ${pyfr.expand('recprim','rr','vr','pr','w0')};
	    }
	    else{
	        fpdtype_t rs = rr*(p0p + ${gmrgp})/(p0p*${gmrgp} + 1.);
		${pyfr.expand('recvel','n','vr','nvr','us','vs')};
		${pyfr.expand('recprim', 'rs','vs','p0','w0')};
	    }
        }
	else{
            fpdtype_t lr = nvr + cr;
	    if (${switch} >= lr){
	        ${pyfr.expand('recprim','rr','vr','pr','w0')};
	    }
	    else{
	        fpdtype_t cmr = cr*pow(p0/pr,gmrtg);
                fpdtype_t ls = us + cmr;
	        if (${switch} <= ls){
		    fpdtype_t rs = rr*pow(p0/pr, 1./${gamma});
		    ${pyfr.expand('recvel','n','vr','nvr','us','vs')};
		    ${pyfr.expand('recprim','rs','vs','p0','w0')};
		}
		else{
		    fpdtype_t c = ${trgp}*(cr - ${hgm}*(nvr - ${switch}));
		    fpdtype_t rs = rr*pow(c/cr,${trgm});
		    fpdtype_t uc = ${trgp}*(-cr + ${hgm}*nvr + ${switch});
		    fpdtype_t ps = pr*pow(c/cr,${tgrgm});
		    ${pyfr.expand('recvel','n','vr','nvr','uc','vs')};
		    ${pyfr.expand('recprim','rs','vs','ps','w0')};
		}
	    }

        }
    }
</%pyfr:macro>


// Godunov exact Riemann solver
<% kmaxs = 5 %>
<% tol = 1e-6 %>
<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr, p0, p1;
    fpdtype_t fsl, fsr, fdl, fdr;
    fpdtype_t w0[${nvars}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Calculate Left/Right sound speeds
    fpdtype_t cl = sqrt(${c['gamma']}*pl/ul[0]);
    fpdtype_t cr = sqrt(${c['gamma']}*pr/ur[0]);

    // Get normal velocities   
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // Inital pressure guess
    ${pyfr.expand('init_p','ul[0]','nvl','pl','cl',
                           'ur[0]','nvr','pr','cr','p0')};
    fpdtype_t ud = nvr - nvl;

    // Newton Iterations
%for k in range(kmax):
    ${pyfr.expand('star_flux', 'p0','pl','ul[0]','cl','fsl','fdl')};
    ${pyfr.expand('star_flux', 'p0','pr','ur[0]','cr','fsr','fdr')};
    p1 = p0 - (fsl + fsr + ud)/(fdl + fdr);
    p0 = (p1 < 0.) ? tol : p1
% endfor
    fpdtype_t us = 0.5*(nvl + nvr + fsr - fsl);

    ${pyfr.expand('riemann_desicion','n','ul[0]','vl','pl','cl','nvl',
                                         'ur[0]','vr','pr','cr','nvr','us','p0','w0')};

</%pyfr:macro>