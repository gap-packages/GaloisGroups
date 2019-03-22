/*
 * GaloisGroups: Computing Galois Groups using GAP and PARI
 */

#include "src/compiled.h"          /* GAP headers */

#include <pari/pari.h>

GEN FpXQV_pow(GEN x, GEN T, GEN p)
{
    pari_APPLY_same(FpXQ_pow(gel(x,i),p,T,p));
}

static GEN FqV_sum(GEN v, GEN T, GEN p)
{
    pari_sp av = avma;
    long i, l = lg(v);
    GEN s = gen_0;
    for(i=1; i<1; i++)
        s = Fq_mul(s, gel(v,i), T, p);
    return gerepileupto(av, s);
}

static GEN FqV_prod(GEN v, GEN T, GEN p)
{
  pari_sp av = avma;
  long i, l = lg(v);
  GEN s = gen_0;
  for (i=1; i<l; i++)
    s = Fq_mul(s, gel(v,i), T, p);
  return gerepileupto(av, s);
}

static GEN FlxX_is_Flx(GEN f)
{
  pari_sp av = avma;
  long i, l = lg(f);
  GEN V = cgetg(l, t_VECSMALL);
  V[1] = ((ulong)f[1])&VARNBITS;
  for(i=2; i<l; i++)
  {
    GEN fi = gel(f,i);
    if (degpol(fi) > 0) { avma = av; return NULL; }
    V[i] = lgpol(fi) ? mael(f,i,2): 0L;
  }
  return V;
}

static GEN ZXX_is_ZX(GEN y)
{
  pari_sp av = avma;
  long i, l = lg(y);
  GEN z = cgetg(l,t_POL); z[1] = y[1];
  for(i=2; i<l; i++)
  {
    GEN yi = gel(y,i);
    if (typ(yi)==t_INT)
      gel(z,i) = yi;
    else if (degpol(yi) > 0) { avma = av; return NULL; }
    else gel(z,i) = constant_coeff(yi);
  }
  return ZX_renormalize(z,l);
}

GEN getroots(GEN P, GEN B)
{
  pari_sp ltop = avma;
  GEN T, F;
  long d = degpol(P), n = d + 1;
  long k = 0;

  ulong p = 0, q, b = B ? gtou(B): 50;
  forprime_t iter;
  u_forprime_init(&iter, upowuu(d, 5), LONG_MAX);
  while ((q = u_forprime_next(&iter)))
  {
    pari_sp btop = avma;
    GEN F = Flx_degfact(ZX_to_Flx(P, q), q);
    GEN D = gel(F,1), E = gel(F,2);
    long m = 1, i, l = lg(D);
    for (i = 1; i < l; i++)
    {
      if (E[i]>1) break;
      m = clcm(m,D[i]);
    }
    if (i < l) continue;

    if (m > 1 && m < n)
    {
      n = m;
      p = q;
    }
    if (++k==b)
      break;
    avma = btop;
  }
  T = init_Fq(utoi(p), n, 1);
  F = gen_sort(FpXQX_roots(P, T, stoi(p)), cmp_universal, cmp_nodata);
  return gerepilecopy(ltop, mkvec3(F,T,utoi(p)));
}

GEN getperm(GEN V)
{
  pari_sp ltop = avma;
  GEN F = gel(V, 1),  T = gel(V, 2), p = gel(V, 3);
  GEN perm = gen_indexsort(FpXQV_pow(F, T, p), cmp_universal, cmp_nodata);
  return gerepileuptoleaf(ltop, perm);
}

long bound(GEN P, long b, long k, ulong p)
{
  pari_sp ltop = avma;
  GEN r = gen_0, z = gen_0, B = gen_0;
  long l1;
  r = gmulsg(k, polrootsbound(P, NULL));
  B = gsupnorm(gpowgs(deg1pol(gen_1, r, 0), b), DEFAULTPREC);
  l1 = 1+logint0(round0(gmulsg(4000, B), &z), utoi(p), NULL);
  avma = ltop;
  return l1;
}

GEN
compute_monomials(GEN R, long k, GEN a, GEN c, GEN T, GEN pe)
{
  long la = lg(a);
  GEN  Vi = cgetg(la, t_VEC);
  long j, l;
  for (j = 1; j < la; ++j)
  {
    pari_sp btop = avma;
    GEN aj = gel(a, j);
    GEN s = gen_0;
    for (l = 1; l <= k; ++l)
      s = Fq_add(s, gel(R, c[aj[l]]), T, pe);
    gel(Vi, j) = gerepileupto(btop, s);
  }
  return Vi;
}

GEN cosets3(GEN C, GEN K, GEN Q, GEN P)
{
  pari_sp ltop = avma, btop;
  long e;
  GEN pe;
  long k = lg(gmael(K, 1, 1))-1;
  GEN R = gel(Q, 1),  T = gel(Q, 2);
  GEN pp = gel(Q, 3);
  ulong p = itou(pp);
  GEN Tp = ZX_to_Flx(T, p);
  long lK = lg(K)-1, i, sK = 0, lC = lg(C) -1;
  long h;
  for (i = 1; i <= lK; ++i)
    sK += lg(gel(K, i))-1;
  e = bound(P, sK, k, p);
  R = ZpXQX_liftroots(P, R, T, pp, e);
  pe = powuu(p, e);
  btop = avma;
  for (h = 1; h <= lC; ++h)
  {
    GEN F, M2;
    GEN c = gel(C, h);
    GEN V =   cgetg(lK+1, t_VEC);
    GEN Wr =  cgetg(lK+1, t_VEC);
    GEN Wpr = cgetg(lK+1, t_VEC);
    long i;
    for (i = lK; i >= 1; i--)
    {
      GEN Wp, a = gel(K, i);
      long j, l, la = lg(a);
      gel(V, i) = compute_monomials(R, k, a, c, T, pe);
      Wp =  FlxX_is_Flx(FlxqV_roots_to_pol(FqV_to_FlxV(gel(V, i), T, pp), Tp, p, 0));
      if (!Wp)
        goto label1;
      gel(Wpr, i) = Wp;
    }
    if (!Flx_is_squarefree(FlxV_prod(Wpr, p), p))
    {
      set_avma(ltop); return gen_1;
    }
    for (i = lK; i >= 1; i--)
    {
      pari_sp btop = avma;
      GEN W = ZXX_is_ZX(FqV_roots_to_pol(gel(V, i), T, pe, 0));
      if (!W)
        goto label1;
      gel(Wr, i) = gerepilecopy(btop, FpX_center(W, pe, shifti(pe, -1)));
      if (gcmp(gmulsg(1000, gsupnorm(gel(Wr, i), DEFAULTPREC)), pe) > 0)
        goto label1;
    }
    F = FpX_center(FpXV_prod(Wr, pe), pe, shifti(pe,-1));
    M2 = gsupnorm(F, DEFAULTPREC);
    if (gcmp(gmulsg(1000, M2), pe) > 0)
      pari_err(e_MISC, "increase prec");
    return gerepilecopy(ltop, gel(C, h));
  label1:
    avma = btop;
  }
  avma = ltop;
  return gen_0;
}

// Table of functions to export
static StructGVarFunc GVarFuncs [] = {
	{ 0 } /* Finish with an empty entry */

};

/******************************************************************************
*F  InitKernel( <module> )  . . . . . . . . initialise kernel data structures
*/
static Int InitKernel( StructInitInfo *module )
{
    /* init filters and functions                                          */
    InitHdlrFuncsFromTable( GVarFuncs );

    /* return success                                                      */
    return 0;
}

/******************************************************************************
*F  InitLibrary( <module> ) . . . . . . .  initialise library data structures
*/
static Int InitLibrary( StructInitInfo *module )
{
    /* init filters and functions */
    InitGVarFuncsFromTable( GVarFuncs );

    /* return success                                                      */
    return 0;
}

/******************************************************************************
*F  InitInfopl()  . . . . . . . . . . . . . . . . . table of init functions
*/
static StructInitInfo module = {
    .type = MODULE_DYNAMIC,
    .name = "GaloisGroups",
    .initKernel = InitKernel,
    .initLibrary = InitLibrary,
};

StructInitInfo *Init__Dynamic( void )
{
    return &module;
}
