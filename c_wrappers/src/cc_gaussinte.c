
// prototype

void f90_gaussinte_w_(int* n, double* x, double* w, double* eps);
double f90_gaussinte_fx_(int* n, double* x, double* a, double* b, double* fx);
double f90_gaussinte_fv_(int* n, double* w, double* c, double* fv);

//

void cc_gaussinte_w(int n, double* x, double* w, double eps) {
  f90_gauusinte_w_(&n, x, w, &eps);
}

double cc_gaussinte_fx(int n, double* x, double a, double b, double* fx) {
  return f90_gaussinte_fx_(&n, x, &a, &b, fx);
}

double cc_gaussinte_fv(int n, double* w, double c, double* fv) {
  return f90_gaussinte_fv_(&n, w, &c, fv);
}
