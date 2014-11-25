#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include "matrix.h"
#include "ilu.h"

using namespace std;

double** rotate(double** A, double* b, int n) {
  for (int i = 0; i < n; i++) {
    double c = A[i][i]/sqrt(A[i][i]*A[i][i] + A[i+1][i]*A[i+1][i]);
    double s = A[i+1][i]/sqrt(A[i][i]*A[i][i] + A[i+1][i]*A[i+1][i]);
    for (int j = i; j < n; j++) {
      double old = A[i][j];
      A[i][j] = old*c + A[i+1][j]*s;
      A[i+1][j] = A[i+1][j]*c - old*s;
    }
    double old = b[i];
    b[i] = old*c;
    b[i+1] = -old*s;
  }
  return A;
}

double* gmres_prec(double** A, double* b, double* x0, int n, int m, double** M) {
  double* r0;
  double betha;
  double** v = new double*[m+1]();
  double** H = new double*[m+1];
  double* f = new double[m+1];
  double** w = new double*[m]();
  double* y = NULL;
  double* xm = new double[n];
  double* tmp1;
  int dim = 0;
  for (int i = 0; i < n; i++) {
    xm[i] = x0[i];
  }
  for (int i = 0; i < m+1; i++) {
    H[i] = new double[m]();
  }
  tmp1 = mult_matrix_to_vector(A, x0, n);
  tmp1 = mult_vector_to_scalar_inplace(tmp1, -1.0, n);
  tmp1 = add_vector_to_vector_inplace(tmp1, b, n);
  r0 = LUSolve(M, tmp1, n);
  f[0] = betha = vector_norm_2(r0, n);
  if (fabs(betha) < 1e-9) {
    goto out;
  }
  v[0] = mult_vector_to_scalar(r0, 1.0/betha, n);
  for (int j = 0; j < m; j++) {
    dim++;
    double* tmp2 = mult_matrix_to_vector(A, v[j], n);
    w[j] = LUSolve(M, tmp2, n);
    delete[] tmp2;
    for (int i = 0; i <= j; i++) {
      double* tmp;
      H[i][j] = mult_vector_to_vector(w[j], v[i], n);
      tmp = mult_vector_to_scalar(v[i], H[i][j], n);
      w[j] = sub_vector_from_vector_inplace(w[j], tmp, n);
      delete[] tmp;
    }
    H[j+1][j] = vector_norm_2(w[j], n);
    if (fabs(H[j+1][j]) < 1e-9) {
      break;
    }
    v[j+1] = mult_vector_to_scalar(w[j], 1.0/H[j+1][j], n);
  }
  H = rotate(H, f, dim);
  y = new double[m];
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
  for (int i = 0; i < m+1; i++) {
    delete[] H[i];
  }
  if (v[0]) for (int i = 0; i < m+1; i++) {
    delete[] v[i];
    if (i < m) delete[] w[i];
  }
  delete[] v;
  delete[] H;
  delete[] f;
  delete[] w;
  if (y) delete[] y;
  delete[] tmp1;

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
  double** M;
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
  A = read_MatrixMarket("tests/watt__1.mtx", n);
  b = new double[n];
  for (int i = 0; i < n; i++) {
    b[i] = (1e-4*rand())/RAND_MAX;
  }
  M = ilu(A, n);
  xm = new double[n]();
  while (diff(A, b, xm, n) > 1e-6) {
    if (cnt >= 300) break;
    cnt++;
    double* x = gmres_prec(A, b, xm, n, 10, M);
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
    delete[] M[i];
  }
  delete[] A;
  delete[] M;
  delete[] b;
  delete[] xm;
  return 0;
}
