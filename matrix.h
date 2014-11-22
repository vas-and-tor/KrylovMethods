#include <cmath>

double mult_vector_to_vector(double* v1, double* v2, int n) {
  double res = 0.0;
  for (int i = 0; i < n; i++) {
    res += v1[i]*v2[i];
  }
  return res;
}

double vector_norm_2(double* v, int n) {
  return sqrt(mult_vector_to_vector(v, v, n));
}

double* mult_matrix_to_vector(double** A, double* v, int n) {
  double* w = new double[n];
  for (int i = 0; i < n; i++) {
    w[i] = mult_vector_to_vector(A[i], v, n);
  }
  return w;
}

double* mult_vector_to_scalar(double* v, double c, int n) {
  double* w = new double[n];
  for (int i = 0; i < n; i++) {
    w[i] = c * v[i];
  }
  return w;
}

double* mult_vector_to_scalar_inplace(double* v, double c, int n) {
  for (int i = 0; i < n; i++) {
    v[i] = c * v[i];
  }
  return v;
}

double* add_vector_to_vector(double* v1, double* v2, int n) {
  double* w = new double[n];
  for (int i = 0; i < n; i++) {
    w[i] = v1[i] + v2[i];
  }
  return w;
}

double* sub_vector_from_vector(double* v1, double* v2, int n) {
  double* w = new double[n];
  for (int i = 0; i < n; i++) {
    w[i] = v1[i] - v2[i];
  }
  return w;
}

double* sub_vector_from_vector_inplace(double* v1, double* v2, int n) {
  for (int i = 0; i < n; i++) {
    v1[i] = v1[i] - v2[i];
  }
  return v1;
}

double* add_vector_to_vector_inplace(double* v1, double* v2, int n) {
  for (int i = 0; i < n; i++) {
    v1[i] = v1[i] + v2[i];
  }
  return v1;
}
