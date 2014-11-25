#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include "matrix.h"
#include "ilu.h"
#include <cmath>

using namespace std;

double diff(double** A, double* b, double* X, int n) {
  double tmp[n];
  mult_matrix_to_vector_inplace(A, X, n, tmp);
  double ans;
  mult_vector_to_scalar_inplace(tmp, -1.0, n);
  add_vector_to_vector_inplace(tmp, b, n);
  ans = vector_norm_2(tmp, n);
  return ans;
}

double * bcgstab_prec(double** A, double* b, double* x0, int n, int & cnt, double** M) {
  double r1[n];
  double r2[n];
  double* tmp = mult_matrix_to_vector(A, x0, n);
  double rho[2];
  double w = 0.0;
  double alpha = 0;
  double v[n] = {};
  double p1[n] = {};
  double* x = new double[n];
  for (int i = 0; i < n; i++) {
    r1[i] = r2[i] = b[i] - tmp[i];
    x[i] = x0[i];
  }
  cnt = 0;
  for (int j = 0; j < 300; j++) {
    cnt++;
    double betha;
    double* p2;
    double s1[n];
    double* s2;
    double t[n];
    rho[j%2] = mult_vector_to_vector(r1, r2, n);
    if (fabs(rho[j%2]) < 1e-15) goto out;
    betha = (j != 0 ? ((rho[j%2]/rho[1-j%2]) * (alpha/w)) : 0);
    for (int i = 0; i < n; i++) {
      p1[i] = r1[i] + betha*(p1[i] - w*v[i]);
    }
    p2 = LUSolve(M, p1, n);
    mult_matrix_to_vector_inplace(A, p2, n, v);
    alpha = rho[j%2]/mult_vector_to_vector(v, r2, n);
    for (int i = 0; i < n; i++) {
      s1[i] = r1[i] - alpha*v[i];
    }
    if (vector_norm_2(s1, n) < 1e-6) {
      for (int i = 0; i < n; i++) {
        x[i] += alpha*p2[i];
      }
      delete[] p2;
      goto out;
    }
    s2 = LUSolve(M, s1, n);
    mult_matrix_to_vector_inplace(A, s2, n, t);
    w = mult_vector_to_vector(s1, t, n)/mult_vector_to_vector(t, t, n);
    for (int i = 0; i < n; i++) {
      x[i] += alpha*p2[i] + w*s2[i];
      r1[i] = s1[i] - w*t[i];
    }
    delete[] p2;
    delete[] s2;
    cerr << j+1 << " Answer: ";
    for (int i = 0; i < n; i++) {
      cerr << x[i] << " ";
    }
    cerr << endl;
    cerr << "diff = " << diff(A, b, x, n) << endl;
    if (vector_norm_2(r1, n) < 1e-6) goto out;
  }

out:

  delete[] tmp;

  return x;
}

int main(void) {
  int n;
  double** A;
  double* b;
  double* xm;
  double* x0;
  double** M;
  int cnt = 0;
  //freopen("input.txt", "rt", stdin);
  //freopen("output.txt", "wt", stdout);
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
  freopen("output.txt", "wt", stdout);
  A = read_MatrixMarket("tests/watt__1.mtx", n);
  b = new double[n];
  for (int i = 0; i < n; i++) {
    b[i] = (1e-4*rand())/RAND_MAX;
  }
  M = ilu(A, n);
  x0 = new double[n];
  for (int i = 0; i < n; i++) {
    x0[i] = 0;
  }
  xm = bcgstab_prec(A, b, x0, n, cnt, M);
  cerr << "Number of iterations: " << cnt << endl;
  cerr << "diff = " << diff(A, b, xm, n) << endl;
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
  delete[] x0;
  delete[] xm;
  return 0;
}

