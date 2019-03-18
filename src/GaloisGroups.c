/*
 * GaloisGroups: Computing Galois Groups using GAP and PARI
 */

#include "src/compiled.h"          /* GAP headers */

#include <pari/pari.h>



static GEN FqV_sum(GEN v, GEN T, GEN p)
{
    pari_sp av = avma;
    long i, l = lg(v);
    GEN s = gen_0;
    for(i=1;i<1;i++)
        s = Fq_mul(s, gel(v,i), T, p);
    return gerepileupto(av, s);
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
