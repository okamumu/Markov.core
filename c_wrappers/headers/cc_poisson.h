#ifndef _CC_POISSON_H
#define _CC_POISSON_H

#ifdef __cplusplus
extern "C" {
#endif

///

int cc_poisson_rightbound(double lambda, double eps);

double cc_poisson_prob(double lambda, int left, int right, double* prob);

#ifdef __cplusplus
}
#endif
#endif
