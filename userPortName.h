/*
 **** ADAMS/Solver include %W% %G%
 */
/*
//
// ***********************************************************
// *****               COPYRIGHT (C) 2004                *****
// *****            BY MSC Software Corporation          *****
// *****              ALL RIGHTS RESERVED                *****
// ***********************************************************
//
//
// Name:
// userPortName.h
//
// Description: C wrapers for FORTRAN routines
//
*/
#ifndef USERPORTNAMES_H
#define USERPORTNAMES_H 1 

#include <string.h>


#ifndef UTILPORTNAMES_H
  typedef int            BOOL;
  typedef double         REAL;
# define UTILPORTNAMES_H
# ifdef FALSE
#   undef FALSE 
# endif
# ifdef TRUE
#   undef TRUE
# endif
# define  TRUE          (BOOL) 1
# define  FALSE         (BOOL) 0
#endif

#ifndef PUBLIC
#ifdef __cplusplus
#define PUBLIC extern "C"
#else
#define PUBLIC
#endif
#endif


/* Handle weird MicroSoft function definitions (FORTRAN to C) */

#if (defined port2pcnt)
#  define BOOL_FUNCTION    PUBLIC BOOL __stdcall
#  define VOID_FUNCTION    PUBLIC void __stdcall
#  define REAL_FUNCTION    PUBLIC REAL __stdcall
#  define INT_FUNCTION     PUBLIC int  __stdcall
#else
#  define BOOL_FUNCTION    PUBLIC BOOL
#  define VOID_FUNCTION    PUBLIC void
#  define REAL_FUNCTION    PUBLIC REAL
#  define INT_FUNCTION     PUBLIC int
#endif


/* Handle function-name conversion to lowercase w/ or w/o trailing underscore */


#if (port2hp700) || (defined port2ibmrs) 

#define   CONSUB     consub
#define   CURSUB     cursub
#define   CFFSUB     cffsub
#define   CNFSUB     cnfsub
#define   COUSUB     cousub
#define   COUXX      couxx 
#define   COUXX2     couxx2
#define   DIFSUB     difsub
#define   DMPSUB     dmpsub
#define   FIESUB     fiesub
#define   GFOSUB     gfosub
#define   GSE_DERIV  gse_deriv
#define   GSE_UPDATE gse_update
#define   GSE_OUTPUT gse_output
#define   GSE_SAMP   gse_samp
#define   GSESUB     gsesub
#define   GSEXU      gsexu
#define   GSEXX      gsexx
#define   GSEYU      gseyu
#define   GSEYX      gseyx
#define   MFOSUB     mfosub
#define   MOTSUB     motsub
#define   RELSUB     relsub
#define   REQSUB     reqsub
#define   SAVSUB     savsub
#define   SENSUB     sensub
#define   SFOSUB     sfosub
#define   SPLINE_READ spline_read
#define   TIRSUB     tirsub
#define   UCOSUB     ucosub
#define   VARSUB     varsub
#define   VFOSUB     vfosub
#define   VTOSUB     vtosub

#endif

#if (defined port2sun) || (port2solaris)  || (defined port2sgi) || (defined port2linux)

#define   CONSUB     consub_
#define   CURSUB     cursub_
#define   CFFSUB     cffsub_
#define   CNFSUB     cnfsub_
#define   COUSUB     cousub_
#define   COUXX      couxx_
#define   COUXX2     couxx2_
#define   DIFSUB     difsub_
#define   DMPSUB     dmpsub_
#define   FIESUB     fiesub_
#define   GFOSUB     gfosub_
#define   GSE_DERIV  gse_deriv_
#define   GSE_UPDATE gse_update_
#define   GSE_OUTPUT gse_output_
#define   GSE_SAMP   gse_samp_
#define   GSESUB     gsesub_
#define   GSEXU      gsexu_
#define   GSEXX      gsexx_
#define   GSEYU      gseyu_
#define   GSEYX      gseyx_
#define   MFOSUB     mfosub_
#define   MOTSUB     motsub_
#define   RELSUB     relsub_
#define   REQSUB     reqsub_
#define   SAVSUB     savsub_
#define   SENSUB     sensub_
#define   SFOSUB     sfosub_
#define   SPLINE_READ spline_read_
#define   TIRSUB     tirsub_
#define   UCOSUB     ucosub_
#define   VARSUB     varsub_
#define   VFOSUB     vfosub_
#define   VTOSUB     vtosub_

#endif


#if (defined port2pcnt)

/*  Function names not modified, they stay upper case as is  */

#endif

/* prototypes for strings */
/* TIRSUB */
#if (defined port2pcnt)
VOID_FUNCTION TIRSUB(int *id, REAL *time, REAL *t0, 
		     REAL *cprop, REAL *tprop, REAL *mprop, 
		     REAL *par, int *npar, char *str, int len, int *nstr, 
		     BOOL *dflag, BOOL *iflag, 
		     REAL *fsae, REAL *tdse, REAL *fprop);

VOID_FUNCTION SPLINE_READ(int *id, char *file, int lfile, 
			  char *block, int lblock, BOOL *errflg);

#elif (defined port2cray)
VOID_FUNCTION TIRSUB(int *id, REAL *time, REAL *t0, 
		     REAL *cprop, REAL *tprop, REAL *mprop, 
		     REAL *par, int *npar, _fcd fcdStr, int *nstr, 
		     BOOL *dflag, BOOL *iflag, 
		     REAL *fsae, REAL *tdse, REAL *fprop);

VOID_FUNCTION SPLINE_READ(int *id, _fcd file, _fcd  block, BOOL *errflg);

#else
VOID_FUNCTION TIRSUB(int *id, REAL *time, REAL *t0, 
		     REAL *cprop, REAL *tprop, REAL *mprop, 
		     REAL *par, int *npar, char *str, int *nstr, 
		     BOOL *dflag, BOOL *iflag, 
		     REAL *fsae, REAL *tdse, REAL *fprop, int len);

VOID_FUNCTION SPLINE_READ(int *id, char *file, char *block, BOOL *errflg, 
			  int lfile, int lblock);

#endif

/* the rest of the prototypes */
  /* CFFSUB */
VOID_FUNCTION CFFSUB(int *id, REAL *time, REAL *par, int *npar, REAL *loci, REAL *locj, REAL *x, 
		     REAL *xdot, REAL *nforce, BOOL *dflag, BOOL *iflag, REAL *result);
  /* CNFSUB */
VOID_FUNCTION CNFSUB(int *id, REAL *time, REAL *par, int *npar, REAL *loci, REAL *ni, REAL *locj, REAL *nj, REAL *gap,
		     REAL *gapdot, REAL *gapdotdot, REAL *area, BOOL *dflag, BOOL *iflag, REAL *result);
  /* COUSUB */
VOID_FUNCTION COUSUB(int *id, REAL *time, REAL *par, int *npar, REAL *disp, int *ndisp, BOOL *iflag, REAL *phi);
  /* COUXX  */
VOID_FUNCTION COUXX (int *id, REAL *time, REAL *par, int *npar, REAL *disp, int *ndisp, BOOL *iflag, REAL *dfda);
  /* COUXX2 */
VOID_FUNCTION COUXX2(int *id, REAL *time, REAL *par, int *npar, REAL *disp, int *ndisp, BOOL *iflag, REAL *dfda2);
  /* CONSUB */
VOID_FUNCTION CONSUB(REAL *par, int *nPar);
/* CURSUB */
VOID_FUNCTION CURSUB(int *id, REAL *par, int *nPar, REAL *alpha, int *iord, BOOL *iflag, REAL *values);
/* DIFSUB */
VOID_FUNCTION DIFSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, REAL *values);
/* DMPSUB */
VOID_FUNCTION DMPSUB(int *id, REAL *time, REAL *par, int *nPar, REAL *freq, int *nmode, REAL *h, REAL *cratios);
/* FIESUB */
VOID_FUNCTION FIESUB(int *id, REAL *time, REAL *par, int *nPar, REAL *disp, REAL *velo, BOOL *dflag, 
		     BOOL *iflag, REAL *field, REAL *dfddis, REAL *dfdvel);
/* GFOSUB */
VOID_FUNCTION GFOSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, REAL *result);
/* GSE_DERIV */
VOID_FUNCTION GSE_DERIV(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, int *ns, REAL *xdot);
/* GSE_UPDATE */
VOID_FUNCTION GSE_UPDATE(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, int *nd, REAL *xdplus1);
/* GSE_OUTPUT */
VOID_FUNCTION GSE_OUTPUT(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, int *no, REAL *y);
/* GSE_SAMP */
VOID_FUNCTION GSE_SAMP(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, REAL *sample_step);
/* GSESUB */
VOID_FUNCTION GSESUB(int *id, REAL *time, REAL *par, int *nPar, 
		     BOOL *dflag, BOOL *iflag, int *nstate, REAL *states, 
		     int *ninput, REAL *input, int *noutpt, REAL *started, 
		     REAL *output);
/* GSEXU */
VOID_FUNCTION GSEXU(int *id, REAL *time, REAL *par, int *nPar, 
		    BOOL *iflag, int *nstate, REAL *states, 
		    int *ninput, REAL *input, int *noutpt, REAL *derivs);
/* GSEXX */
VOID_FUNCTION GSEXX(int *id, REAL *time, REAL *par, int *nPar, 
		    BOOL *iflag, int *nstate, REAL *states, 
		    int *ninput, REAL *input, int *noutpt, REAL *derivs);
/* GSEYU */
VOID_FUNCTION GSEYU(int *id, REAL *time, REAL *par, int *nPar, 
		    BOOL *iflag, int *nstate, REAL *states, 
		    int *ninput, REAL *input, int *noutpt, REAL *derivs);
/* GSEYX */
VOID_FUNCTION GSEYX(int *id, REAL *time, REAL *par, int *nPar, 
		    BOOL *iflag, int *nstate, REAL *states, 
		    int *ninput, REAL *input, int *noutpt, REAL *derivs);
/* MFOSUB */
VOID_FUNCTION MFOSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, 
                     REAL *modloads, int *nmodes, int *ncases, REAL *scale, int *icase, REAL *loadvec);
/* MOTSUB */
VOID_FUNCTION MOTSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, REAL *value);
/* RELSUB */
VOID_FUNCTION RELSUB(int *iunit, BOOL *errflg);
/* REQSUB */
VOID_FUNCTION REQSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *iflag, REAL *result);
/* SAVSUB */
VOID_FUNCTION SAVSUB(int *iunit, BOOL *errflg);
/* SENSUB */
VOID_FUNCTION SENSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *iflag, REAL *value);
/* SEVSUB */
VOID_FUNCTION SEVSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *iflag, REAL *result);
/* SFOSUB */
VOID_FUNCTION SFOSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, REAL *value);

/* UCOSUB */
VOID_FUNCTION UCOSUB(int *id, REAL *time, REAL *q, REAL *par, int *npar, 
		     int *idrsel, BOOL *iflag, 
		     REAL *scalar, REAL *array, REAL *xmatrx);
/* VARSUB */
VOID_FUNCTION VARSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, REAL *value);
/* VFOSUB */
VOID_FUNCTION VFOSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, REAL *result);
/* VTOSUB  */
VOID_FUNCTION VTOSUB(int *id, REAL *time, REAL *par, int *nPar, BOOL *dflag, BOOL *iflag, REAL *result);

#endif
