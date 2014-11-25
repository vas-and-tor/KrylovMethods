#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "matrix.h"
#include <cmath>
#include <cstdio>

using namespace std;

double diff(double** A, double* b, double* X, int n) {
  double* tmp = mult_matrix_to_vector(A, X, n);
  double ans;
  tmp = mult_vector_to_scalar_inplace(tmp, -1.0, n);
  tmp = add_vector_to_vector_inplace(tmp, b, n);
  ans = vector_norm_2(tmp, n);
  delete[] tmp;
  return ans;
}

double * bcg(double** A, double* b, double* x0, int n, int & cnt) {
  double r1[n];
  double r2[n];
  double** AT = new double*[n];
  double* tmp = mult_matrix_to_vector(A, x0, n);
  double rho[2];
  double p1[n] = {};
  double p2[n] = {};
  double* x = new double[n];
  for (int i = 0; i < n; i++) {
    r1[i] = r2[i] = b[i] - tmp[i];
    AT[i] = new double[n];
    x[i] = x0[i];
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      AT[i][j] = A[j][i];
    }
  }
  cnt = 0;
  for (int j = 0; j < 300; j++) {
    cnt++;
    double q1[n];
    double q2[n];
    double betha;
    double alpha;
    rho[j%2] = mult_vector_to_vector(r2, r1, n);
    if (fabs(rho[j%2]) < 1e-15) {
      goto out;
    }
    betha = (j != 0 ? rho[j%2]/rho[1-j%2] : 0);
    for (int i = 0; i < n; i++) {
      p1[i] = r1[i] + betha*p1[i];
      p2[i] = r2[i] + betha*p2[i];
    }
    mult_matrix_to_vector_inplace(A , p1, n, q1);
    mult_matrix_to_vector_inplace(AT, p2, n, q2);
    alpha = rho[j%2]/mult_vector_to_vector(q1, p2, n);
    for (int i = 0; i < n; i++) {
      x[i] += alpha*p1[i];
      r1[i] -= alpha*q1[i];
      r2[i] -= alpha*q2[i];
    }
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
  for (int i = 0; i < n; i++) {
    delete[] AT[i];
  }
  delete[] AT;

  return x;
}

int main(void) {
  int n;
  double** A;
  double* b;
  double* xm;
  double* x0;
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
    b[i] = (1e4*rand())/RAND_MAX;
  }
  x0 = new double[n];
  for (int i = 0; i < n; i++) {
    x0[i] = 0;
  }
  xm = bcg(A, b, x0, n, cnt);
  cerr << "Number of iterations: " << cnt << endl;
  cerr << "diff = " << diff(A, b, xm, n) << endl;
  for (int i = 0; i < n; i++) {
    cout << fixed << setprecision(5) << xm[i] << endl;
  }
  for (int i = 0; i < n; i++) {
    delete[] A[i];
  }
  delete[] A;
  delete[] b;
  delete[] x0;
  delete[] xm;
  return 0;
}

