/*
 * GaloisGroups: Computing Galois Groups using GAP and PARI
 */

#include "src/compiled.h"          /* GAP headers */

#include <pari/pari.h>

GEN FpXQV_pow(GEN x, GEN T, GEN p)
{
    pari_APPLY_same(FpXQ_pow(gel(x,i),p,T,p));
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

static GEN
Fq_compute_monomials(GEN R, long k, GEN a, GEN c, GEN T, GEN p)
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
      s = Fq_add(s, gel(R, c[aj[l]]), T, p);
    gel(Vi, j) = gerepileupto(btop, s);
  }
  return Vi;
}

static GEN
Flxq_compute_monomials(GEN R, long k, GEN a, GEN c, GEN T, ulong p)
{
  long la = lg(a);
  GEN  g0 = pol0_Flx(get_Flx_var(T));
  GEN  Vi = cgetg(la, t_VEC);
  long j, l;
  for (j = 1; j < la; ++j)
  {
    pari_sp btop = avma;
    GEN aj = gel(a, j);
    GEN s = g0;
    for (l = 1; l <= k; ++l)
      s = Flx_add(s, gel(R, c[aj[l]]), p);
    gel(Vi, j) = gerepileupto(btop, s);
  }
  return Vi;
}

static GEN
FlxqV_prod(GEN v, GEN T, ulong p)
{
  pari_sp av = avma;
  long i, l = lg(v);
  GEN s = pol1_Flx(get_Flx_var(T));
  for (i = 1; i < l; i++)
    s = Flxq_mul(s, gel(v,i), T, p);
  return gerepileupto(av, s);
}

GEN
cosets_squarefree(GEN C, GEN K, GEN Q, GEN P)
{
  pari_sp ltop = avma, btop;
  long k = lg(gmael(K, 1, 1))-1;
  GEN R = gel(Q, 1),  T = gel(Q, 2);
  GEN pp = gel(Q, 3);
  ulong p = itou(pp);
  GEN Tp = ZX_to_Flx(T, p), Rp = ZXV_to_FlxV(R, p);
  long lK = lg(K)-1, lC = lg(C) -1;
  long h;
  btop = avma;
  for (h = 1; h <= lC; ++h)
  {
    GEN c = gel(C, h);
    GEN Wpr = cgetg(lK+1, t_VEC);
    GEN Wps;
    long i;
    for (i = 1; i <= lK; i++)
    {
      GEN V = Flxq_compute_monomials(Rp, k, gel(K, i), c, Tp, p);
      gel(Wpr, i) = FlxqV_prod(V, Tp, p);
    }
    Wps = gtoset(Wpr);
    if (lg(Wps)!=lg(Wpr))
    {
      avma = ltop; return gen_1;
    }
    avma = btop;
  }
  avma = ltop;
  return gen_0;
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
    long i;
    for (i = lK; i >= 1; i--)
      gel(V, i) = Fq_compute_monomials(R, k, gel(K, i), c, T, pe);
    for (i = lK; i >= 1; i--)
    {
      GEN W = ZXX_is_ZX(FqV_roots_to_pol(gel(V, i), T, pe, 0));
      if (!W)
        goto label1;
      gel(Wr, i) = FpX_center(W, pe, shifti(pe, -1));
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
