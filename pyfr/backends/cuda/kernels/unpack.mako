# -*- coding: utf-8 -*-
<%inherit file='base'/>

__global__ void
unpack_view(int n, int nrv, int ncv,
          const fpdtype_t* __restrict__ pmat,
          const int* __restrict__ vix,
          const int* __restrict__ vrstri,
          fpdtype_t* __restrict__ v)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;

    if (i < n && ncv == 1)
        v[vix[i]] = pmat[i];
    else if (i < n && nrv == 1)
        for (int c = 0; c < ncv; ++c)
            v[vix[i] + SOA_SZ*c] = pmat[c*n + i];
    else if (i < n)
        for (int r = 0; r < nrv; ++r)
            for (int c = 0; c < ncv; ++c)
                v[vix[i] + vrstri[i]*r + SOA_SZ*c] = pmat[(r*ncv + c)*n + i];
}
