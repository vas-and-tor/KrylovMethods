#include <cmath>
#include <cstdio>
#include <cstdlib>

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

double** read_MatrixMarket(const char* filename, int & n) {
  double** A;
  FILE* f = fopen(filename, "rt");
  ssize_t read;
  size_t len;
  char* line = NULL;
  n = -1;
  if (f == NULL) {
    fprintf(stderr, "File '%s' not found.", filename);
    exit(1);
  }
  while ((read = getline(&line, &len, f)) != -1) {
    if (line[0] != '%') {
      if (n != -1) {
        int i,j;
        double x;
        sscanf(line, "%d %d %lf", &i, &j, &x);
        A[i-1][j-1] = x;
      } else {
        int nn, mm, ll;
        sscanf(line, "%d %d %d", &nn, &mm, &ll);
        n = (nn >= mm ? nn : mm);
        A = new double*[n];
        for (int i = 0; i < n; i++) {
          A[i] = new double[n];
        }
      }
    }
  }
  return A;
}

double* LUSolve(double** A, double* b, int n) {
  double* tmp = new double[n];
  double* x = new double[n];
  for (int i = 0; i < n; i++) {
    tmp[i] = b[i];
    for (int j = 0; j < i; j++) {
      tmp[i] -= A[i][j]*tmp[j];
    }
  }
  for (int i = n-1; i >= 0; i--) {
    x[i] = tmp[i]/A[i][i];
    for (int j = i+1; j < n; j++) {
      x[i] -= A[i][j]*x[j]/A[i][i];
    }
  }
  delete[] tmp;
  return x;
}
