#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include "matrix.h"

using namespace std;

double * fom(double** A, double* b, double* x0, int n, int m) {
  double* r0;
  double betha;
  double** v = new double*[m];
  double** H = new double*[m];
  double* f = NULL;
  double** w = new double*[m];
  double* y = NULL;
  double* xm = new double[n];
  int dim = 1;
  for (int i = 0; i < n; i++) {
    xm[i] = x0[i];
  }
  for (int i = 0; i < m; i++) {
    H[i] = new double[m]();
  }
  r0 = mult_matrix_to_vector(A, x0, n);
  r0 = mult_vector_to_scalar_inplace(r0, -1.0, n);
  r0 = add_vector_to_vector_inplace(r0, b, n);
  betha = vector_norm_2(r0, n);
  if (fabs(betha) < 1e-9) goto out;
  v[0] = mult_vector_to_scalar(r0, 1.0/betha, n);
  for (int j = 0; j < m; j++) {
    w[j] = mult_matrix_to_vector(A, v[j], n);
    for (int i = 0; i <= j; i++) {
      double* tmp;
      H[i][j] = mult_vector_to_vector(w[j], v[i], n);
      tmp = mult_vector_to_scalar(v[i], H[i][j], n);
      w[j] = sub_vector_from_vector_inplace(w[j], tmp, n);
      delete[] tmp;
    }
    if (j+1 < m) {
      H[j+1][j] = vector_norm_2(w[j], n);
      if (fabs(H[j+1][j]) < 1e-9) {
        break;
      }
      v[j+1] = mult_vector_to_scalar(w[j], 1.0/H[j+1][j], n);
      dim++;
    }
  }
  f = new double[dim]();
  f[0] = betha;
  for (int i = 1; i < dim; i++) {
    for (int j = i; j < dim; j++) {
      H[i][j] -= H[i-1][j]*H[i][i-1]/H[i-1][i-1];
    }
    f[i] -= f[i-1]*H[i][i-1]/H[i-1][i-1];
    H[i][i-1] = 0.0;
  }
  y = new double[dim];
  for (int i = dim-1; i >= 0; i--) {
    y[i] = f[i]/H[i][i];
    for (int j = i+1; j < dim; j++) {
      y[i] -= y[j]*H[i][j]/H[i][i];
    }
  }
  for (int i = 0; i < dim; i++) {
    double* tmp = mult_vector_to_scalar(v[i], y[i], n);
    xm = add_vector_to_vector_inplace(xm, tmp, n);
    delete[] tmp;
  }

out:

  delete[] r0;
  for (int i = 0; i < m; i++) {
    delete[] H[i];
  }
  for (int i = 0; i < dim; i++) {
    delete[] v[i];
    delete[] w[i];
  }
  delete[] v;
  delete[] H;
  if (f) delete[] f;
  delete[] w;
  if (y) delete[] y;

  return xm;
}

double diff(double** A, double* b, double* X, int n) {
  double* tmp = mult_matrix_to_vector(A, X, n);
  double ans;
  tmp = mult_vector_to_scalar_inplace(tmp, -1.0, n);
  tmp = add_vector_to_vector_inplace(tmp, b, n);
  ans = vector_norm_2(tmp, n);
  delete[] tmp;
  return ans;
}

int main(void) {
  int n;
  double** A;
  double* b;
  double* xm;
  int cnt = 0;
  freopen("input.txt", "rt", stdin);
  freopen("output.txt", "wt", stdout);
  //cin >> n;
  //A = new double*[n];
  //for (int i = 0; i < n; i++) {
    //A[i] = new double[n];
  //}
  //b = new double[n];
  //for (int i = 0; i < n; i++) {
    //for (int j = 0; j < n; j++) {
      //cin >> A[i][j];
    //}
  //}
  //for (int i = 0; i < n; i++) {
    //cin >> b[i];
  //}
  A = read_MatrixMarket("cage10.mtx", n);
  b = new double[n];
  for (int i = 0; i < n; i++) {
    b[i] = (1e-4*rand())/RAND_MAX;
  }
  xm = new double[n]();
  while (diff(A, b, xm, n) > 1e-6) {
    if (cnt >= 300) break;
    cnt++;
    double* x = fom(A, b, xm, n, 10);
    for (int i = 0; i < n; i++) {
      xm[i] = x[i];
    }
    cerr << "Answer: ";
    for (int i = 0; i < n; i++) {
      cerr << xm[i] << " ";
    }
    cerr << endl;
    cerr << "diff = " << diff(A, b, xm, n) << endl;
    delete[] x;
  }
  cerr << "Number of iterations: " << cnt << endl;
  for (int i = 0; i < n; i++) {
    cout << fixed << setprecision(5) << xm[i] << endl;
  }
  for (int i = 0; i < n; i++) {
    delete[] A[i];
  }
  delete[] A;
  delete[] b;
  delete[] xm;
  return 0;
}

