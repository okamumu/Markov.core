///

void f90_poisson_rightbound_(double* lambda, double* eps, int* right);

void f90_poisson_prob_(double* lambda, int* left, int* right, double* prob, double* weight);

//

int cc_poisson_rightbound(double lambda, double eps) {
  int right;
  f90_poisson_rightbound_(&lambda, &eps, &right);
  return right;
}

double cc_poisson_prob(double lambda, int left, int right, double* prob) {
  double weight;
  f90_poisson_prob_(&lambda, &left, &right, prob, &weight);
  return weight;
}
