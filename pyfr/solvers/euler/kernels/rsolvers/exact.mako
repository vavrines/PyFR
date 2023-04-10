<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<% gmrtg = (c['gamma'] - 1)/(2*c['gamma']) %>
<% gprtg = (c['gamma'] + 1)/(2*c['gamma']) %>
<% tgrgm = (2*c['gamma'])/(c['gamma'] - 1) %>
<% grgm = c['gamma']/(c['gamma'] - 1) %>
<% trgm = 2/(c['gamma'] - 1) %>
<% trgp = 2/(c['gamma'] + 1) %>
<% gmrgp = (c['gamma'] - 1)/(c['gamma'] + 1) %>
<% hgm = 0.5*(c['gamma'] - 1) %>
<% gamma = c['gamma'] %>
<% p_min = 1e-6 %>

<%pyfr:macro name='init_p' params='rl, vl, pl, cl, rr, vr, pr, cr, p0'>
    // Initial guess for pressure
    fpdtype_t bpv = 0.5*(pl + pr) + 0.125*(vl[0] - vr[0])*(rl + rr)*(cl + cr);
    bpv = max(0.0, bpv);
    fpdtype_t pmin = min(pl, pr);
    fpdtype_t pmax = max(pl, pr);

    if (pmax <= 2*pmin && pmin <= bpv && bpv <= pmax)
    {
        p0 = bpv;
    }
    else if (bpv < pmin)
    {   // Two-rarefaction case
        fpdtype_t pre = pow(pl / pr, ${gmrtg});
        fpdtype_t um  = (pre*vl[0]*cr + vr[0]*cl + ${trgm}*(pre - 1)*cl*cr) /
                        (pre*cr + cl);
        fpdtype_t ptl = 1 - ${hgm}*(um - vl[0]) / cl;
        fpdtype_t ptr = 1 + ${hgm}*(um - vr[0]) / cr;

        p0 = 0.5*(pl*pow(ptl, ${tgrgm}) + pr*pow(ptr, ${tgrgm}));
    }
    else
    {   // Two-shock case
        fpdtype_t gl = sqrt(${trgp} / (rl*(${gmrgp}*pl + bpv)));
        fpdtype_t gr = sqrt(${trgp} / (rr*(${gmrgp}*pr + bpv)));
        p0 = (gl*pl + gr*pr - vr[0] + vl[0]) / (gl + gr);
    }
</%pyfr:macro>

<%pyfr:macro name='star_flux' params='p, ps, rs, cs, f, fd'>
    // Star Flux, assuming covolume = 0. See Toro 2009 Eq.(4.86 - 4.87)
    if (p <= ps)
    {
       fpdtype_t pr = p / ps;
       f = ${trgm}*cs*(pow(pr, ${gmrtg}) - 1);
       fd = pow(pr, ${-gprtg}) / (rs*cs);
    }
    else
    {
       fpdtype_t pbs_inv = 1 / (p + ${gmrgp}*ps);
       fpdtype_t sapb = sqrt(${trgp} * pbs_inv / rs);
       f = (p - ps)*sapb;
       fd = sapb - 0.5*f*pbs_inv;
    }
</%pyfr:macro>

<%pyfr:macro name='riemann_decision' params='rl, vl, pl, cl, rr, vr, pr, cr, us, p0, w0'>
    // Exact solve solution decision tree
    if (0 <= us)
    {
% for i in range(1, ndims):
        w0[${i + 1}] = vl[${i}];
% endfor
        if (p0 <= pl)
        {
            if (vl[0] - cl >= 0)
            {
                w0[0] = rl;
                w0[1] = vl[0];
                w0[${nvars - 1}] = pl;
            }
            else
            {
                fpdtype_t cml = cl*pow(p0 / pl, ${gmrtg});
                if (us - cml < 0)
                {
                    w0[0] = rl*pow(p0 / pl, ${1 / gamma});
                    w0[1] = us;
                    w0[${nvars - 1}] = p0;
                }
                else
                {
                    fpdtype_t c = ${trgp} + ${gmrgp}*vl[0] / cl;
                    w0[0] = rl*pow(c, ${trgm});
                    w0[1] = ${trgp}*cl + ${gmrgp}*vl[0];
                    w0[${nvars - 1}] = pl*pow(c, ${tgrgm});
                }
            }
        }
        else
        {
            fpdtype_t p0p = p0 / pl;
            fpdtype_t sl = vl[0] - cl*sqrt(${gprtg}*p0p + ${gmrtg});
            if (0 <= sl)
            {
                w0[0] = rl;
                w0[1] = vl[0];
                w0[${nvars - 1}] = pl;
            }
            else
            {
                w0[0] = rl*(p0p + ${gmrgp}) / (p0p*${gmrgp} + 1);
                w0[1] = us;
                w0[${nvars - 1}] = p0;
            }
        }
    }
    else
    {
% for i in range(1, ndims):
        w0[${i + 1}] = vr[${i}];
% endfor
        if (p0 > pr)
        {
            fpdtype_t p0p = p0 / pr;
            fpdtype_t sr = vr[0] + cr*sqrt(${gprtg}*p0p + ${gmrtg});
            if (sr <= 0)
            {
                w0[0] = rr;
                w0[1] = vr[0];
                w0[${nvars - 1}] = pr;
            }
            else
            {
                w0[0] = rr*(p0p + ${gmrgp}) / (p0p*${gmrgp} + 1);
                w0[1] = us;
                w0[${nvars - 1}] = p0;
            }
        }
        else
        {
            if (vr[0] + cr <= 0)
            {
                w0[0] = rr;
                w0[1] = vr[0];
                w0[${nvars - 1}] = pr;
            }
            else
            {
                fpdtype_t p0p = p0 / pr;
                fpdtype_t cmr = cr*pow(p0p, ${gmrtg});
                if (us + cmr >= 0)
                {
                    w0[0] = rr*pow(p0p, ${1 / gamma});
                    w0[1] = us;
                    w0[${nvars - 1}] = p0;
                }
                else
                {
                    fpdtype_t c = ${trgp} - ${gmrgp}*vr[0] / cr;
                    w0[0] = rr*pow(c, ${trgm});
                    w0[1] = ${gmrgp}*vr[0] - ${trgp}*cr;
                    w0[${nvars - 1}] = pr*pow(c, ${tgrgm});
                }
            }
        }
    }
</%pyfr:macro>

<%pyfr:macro name='primitives' params='s, p, v, c'>
    fpdtype_t invrho = 1/s[0], E = s[${nvars - 1}];

    // Compute the velocities
% for i in range(ndims):
    v[${i}] = invrho*s[${i + 1}];
% endfor

    // Compute the pressure
    p = ${gamma - 1}*(E - 0.5*invrho*${pyfr.dot('s[{i}]', i=(1, ndims + 1))});

    // Compute the local sound speed
    c = sqrt(${gamma}*p*invrho);
</%pyfr:macro>

<% kmax = 3 %>

<%pyfr:macro name='rsolve_1d' params='ul, ur, nf'>
    fpdtype_t vl[${ndims}], pl, cl, fsl, fdl;
    fpdtype_t vr[${ndims}], pr, cr, fsr, fdr;
    fpdtype_t p0, p1;
    fpdtype_t w0[${nvars}];

    // Compute the left/right srimitives
    ${pyfr.expand('primitives', 'ul', 'pl', 'vl', 'cl')};
    ${pyfr.expand('primitives', 'ur', 'pr', 'vr', 'cr')};

    // Inital pressure guess
    fpdtype_t rl = ul[0];
    fpdtype_t rr = ur[0];
    ${pyfr.expand('init_p', 'rl', 'vl', 'pl', 'cl', 'rr', 'vr', 'pr', 'cr',
                  'p0')};

    // Newton Iterations
% for k in range(kmax):
    // Iteration ${k}
    ${pyfr.expand('star_flux', 'p0', 'pl', 'rl', 'cl', 'fsl', 'fdl')};
    ${pyfr.expand('star_flux', 'p0', 'pr', 'rr', 'cr', 'fsr', 'fdr')};
    p1 = p0 - (fsl + fsr + vr[0] - vl[0]) / (fdl + fdr);
    p0 = (p1 < 0) ? ${p_min} : p1;
% endfor
    fpdtype_t us = 0.5*(vl[0] + vr[0] + fsr - fsl);

    // Riemann solve decision tree
    ${pyfr.expand('riemann_decision', 'rl', 'vl', 'pl', 'cl', 'rr', 'vr', 'pr',
                  'cr', 'us', 'p0', 'w0')};

    // Primitive to 1D flux
    nf[0] = w0[0]*w0[1];
    nf[1] = nf[0]*w0[1] + w0[${nvars - 1}];
% for i in range(2, ndims+1):
    nf[${i}] = nf[0]*w0[${i}];
% endfor
    nf[${nvars - 1}] = (${grgm}*w0[${nvars - 1}] + 0.5*w0[0]*(${pyfr.dot('w0[{i}]', i=(1, ndims + 1))}))*w0[1];

</%pyfr:macro>

<%include file='pyfr.solvers.euler.kernels.rsolvers.rsolve1d'/>
