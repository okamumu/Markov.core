#ifndef _CC_GAUSSINTE_H
#define _CC_GAUSSINTE_H

#ifdef __cplusplus
extern "C" {
#endif

void cc_gaussinte_w(int n, double* x, double* w, double eps);
double cc_gaussinte_fx(int n, double* x, double a, double b, double* fx);
double cc_gaussinte_fv(int n, double* w, double c, double* fv);

#ifdef __cplusplus
}
#endif
#endif
