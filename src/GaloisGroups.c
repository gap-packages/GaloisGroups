/*
 * GaloisGroups: Computing Galois Groups using GAP and PARI
 */

#include "src/compiled.h"          /* GAP headers */

#include <pari/pari.h>

static GEN
FpXQV_pow(GEN x, GEN T, GEN p)
{
    pari_APPLY_same(FpXQ_pow(gel(x,i),p,T,p));
}

GEN
ZXX_is_ZX(GEN y)
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

GEN
getroots(GEN P, GEN B)
{
  pari_sp ltop = avma;
  GEN T, F;
  long d = degpol(P), n = d + 1;
  long k = 0;

  ulong p = 0, q, b = B ? gtou(B): 50;
  forprime_t iter;
  u_forprime_init(&iter, upowuu(d, d), LONG_MAX);
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
  F = FpXQX_roots(P, T, stoi(p));
  return gerepilecopy(ltop, mkvec3(F,T,utoi(p)));
}

GEN
getperm(GEN V)
{
  pari_sp ltop = avma;
  GEN F = gel(V, 1),  T = gel(V, 2), p = gel(V, 3);
  GEN p1 = gen_indexsort(F, cmp_universal, cmp_nodata);
  GEN p2 = gen_indexsort(FpXQV_pow(F, T, p), cmp_universal, cmp_nodata);
  return gerepileuptoleaf(ltop, perm_mul(p2, perm_inv(p1)));
}

static GEN
Fq_compute_monomial(GEN R, long k, GEN K, GEN sig, GEN T, GEN p)
{
  pari_sp av = avma;
  long j, i, l = lg(K);
  GEN  V = gen_0;
  for (j = 1; j < l; ++j)
  {
    GEN Kj = gel(K, j);
    GEN s = gen_1;
    for (i = 1; i <= k; ++i)
      s = Fq_mul(s, gel(R, sig[Kj[i]]), T, p);
    V = Fq_add(V, s, T, p);
    if (gc_needed(av,1))
      V = gerepileupto(av, V);
  }
  return gerepileupto(av, V);
}

static GEN
Fq_compute_monomials(GEN R, long k, GEN K, GEN sig, GEN T, GEN p)
{
  long j, l = lg(sig);
  GEN  V = cgetg(l, t_VEC);
  for (j = 1; j < l; ++j)
    gel(V, j) = Fq_compute_monomial(R, k, K, gel(sig, j), T, p);
  return V;
}

static GEN
Flxq_compute_monomial(GEN R, long k, GEN K, GEN sig, GEN T, ulong p)
{
  pari_sp av = avma;
  long j, i, l = lg(K), v = get_Flx_var(T);
  GEN  V = pol0_Flx(v);
  for (j = 1; j < l; ++j)
  {
    GEN Kj = gel(K, j);
    GEN s = pol1_Flx(v);
    for (i = 1; i <= k; ++i)
      s = Flxq_mul(s, gel(R, sig[Kj[i]]), T, p);
    V = Flxq_add(V, s, T, p);
    if (gc_needed(av,1))
      V = gerepileupto(av, V);
  }
  return gerepileupto(av, V);
}

static GEN
Flxq_compute_monomials(GEN R, long k, GEN K, GEN sig, GEN T, ulong p)
{
  long j, l = lg(sig);
  GEN  V = cgetg(l, t_VEC);
  for (j = 1; j < l; ++j)
    gel(V, j) = Flxq_compute_monomial(R, k, K, gel(sig, j), T, p);
  return V;
}

GEN
cosets_squarefree(GEN C, GEN K, GEN Q, GEN P)
{
  pari_sp ltop = avma;
  long k = lg(gel(K, 1))-1;
  GEN R = gel(Q, 1),  T = gel(Q, 2), pp = gel(Q, 3);
  ulong p = itou(pp);
  GEN Tp = ZX_to_Flx(T, p), Rp = ZXV_to_FlxV(R, p);
  GEN V = Flxq_compute_monomials(Rp, k, K, C, Tp, p);
  if (lg(gtoset(V))!=lg(V))
  {
    avma = ltop; return gen_0;
  }
  avma = ltop; return gen_1;
}

long
bound4(GEN P, long b, long k, long c, ulong p)
{
  pari_sp ltop = avma;
  GEN r = gen_0, z = gen_0, B = gen_0;
  long l1;
  r = gmulgs(gpowgs(polrootsbound(P, NULL), k), b);
  B = gsupnorm(gpowgs(deg1pol(gen_1, r, 0), c), DEFAULTPREC);
  l1 = 1+logint0(round0(gmulsg(1L<<20, B), &z), utoi(p), NULL);
  avma = ltop;
  return l1;
}

GEN
cosets4(GEN C, GEN K, GEN Q, GEN P)
{
  pari_sp av = avma;
  pari_timer ti;
  GEN R = gel(Q, 1),  T = gel(Q, 2);
  GEN pp = gel(Q, 3);
  ulong p = itou(pp);
  GEN pe, pe2,pe20;
  long lK = lg(K)-1, i, lV;
  long k = lg(gel(K, 1))-1;
  GEN V;
  if (DEBUGLEVEL) timer_start(&ti);
  long e = bound4(P, lK, k,lg(C)-1, p);
  R = ZpXQX_liftroots(P, R, T, pp, e);
  if (DEBUGLEVEL) timer_printf(&ti, "lift to %ld",e);
  pe = powuu(p, e); pe2 = shifti(pe,-1); pe20 = shifti(pe,-20);
  V = Fq_compute_monomials(R, k, K, C, T, pe);
  if (DEBUGLEVEL) timer_printf(&ti, "monomials: %ld×%ld×%ld",lK, k,lg(C)-1);
//  GEN W = ZXX_is_ZX(FqV_roots_to_pol(V, T, pe, 0));
//  timer_printf(&ti, "roots_to_pol: %ld",degpol(W));
//  W = FpX_center(W, pe, pe2);
//  if (!ZX_is_squarefree(W)) { avma = av; return gen_1; }
//  timer_printf(&ti, "squarefree");
  lV = lg(V);
  for (i=1; i<lV; i++)
  {
    GEN Vi = gel(V, i);
    if (typ(Vi)==t_POL)
    {
      if (degpol(Vi) > 0) continue;
      Vi = constant_coeff(Vi);
    }
    Vi =  Fp_center(Vi, pe, pe2);
    if (abscmpii(Vi,pe20)>=0) continue;
    avma = av; return gel(C,i);
   // if (signe(poleval(W, Vi))==0) { avma = av; return gel(C,i); }
  }
  if (DEBUGLEVEL) timer_printf(&ti, "poleval: %ld", lV-1);
  avma = av; return gen_0;
}

static GEN
Flx_Flxq_veceval(GEN S, GEN R, GEN T, long p)
{
  long j, l = lg(R);
  GEN  V = cgetg(l, t_VEC);
  for (j = 1; j < l; ++j)
    gel(V, j) = Flx_Flxq_eval(S, gel(R,j), T, p);
  return V;
}

GEN
Tschirnhausen(GEN P, GEN Q)
{
  pari_sp av = avma, av2;
  long v = varn(P);
  GEN R = gel(Q, 1),  T = gel(Q, 2), pp = gel(Q, 3);
  ulong p = itou(pp);
  GEN Tp = ZX_to_Flx(T, p), Rp = ZXV_to_FlxV(R, p);
  GEN S, V, U;
  av2 = avma;
  do
  {
    avma = av2;
    S = random_Flx(4, v, p);
    V = Flx_Flxq_veceval(S, Rp, Tp, p);
  } while (lg(gtoset(V))!=lg(V));
  S = FpX_center(Flx_to_ZX(S), pp, shifti(pp,1));
  U = ZXQ_charpoly(S,P,varn(P));
  V = FlxV_to_ZXV(V);
  Q = mkvec3(V, T, pp);
  return gerepilecopy(av, mkvec2(U, Q));
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
