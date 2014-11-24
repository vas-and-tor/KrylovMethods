#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "matrix.h"
#include <cstdio>

using namespace std;

double * bcgstab(double** A, double* b, double* x0, int n, int & cnt) {
  double* r1[2] = { NULL, new double[n] };
  double* r2 = new double[n];
  double* p = new double[n];
  double* x = new double[n];
  r1[0] = mult_matrix_to_vector(A, x0, n);
  for (int i = 0; i < n; i++) {
    p[i] = r2[i] = r1[0][i] = b[i] - r1[0][i];
    x[i] = x0[i];
  }
  if (vector_norm_2(r1[0], n) < 1e-9) {
    goto out;
  }
  for (int j = 0; j < 300; j++) {
    double* tmp1 = mult_matrix_to_vector(A , p, n);
    double* tmp2;
    double alpha = mult_vector_to_vector(r1[j%2], r2, n)/mult_vector_to_vector(tmp1, r2, n);
    double betha;
    double s[n];
    double w;
    for (int i = 0; i < n; i++) {
      s[i] = r1[j%2][i] - alpha*tmp1[i];
    }
    tmp2 = mult_matrix_to_vector(A, s, n);
    w = mult_vector_to_vector(tmp1, s, n)/mult_vector_to_vector(tmp1, tmp1, n);
    for (int i = 0; i < n; i++) {
      x[i] += alpha*p[i] + w*s[i];
      r1[1-j%2][i] = s[i] - w*tmp2[i];
    }
    betha = mult_vector_to_vector(r1[1-j%2], r2, n)/mult_vector_to_vector(r1[j%2], r2, n) * alpha/w;
    for (int i = 0; i < n; i++) {
      p[i] = r1[1-j%2][i] + betha*(p[i] - w*tmp1[i]);
    }
    cerr << j+1 << " Answer: ";
    for (int i = 0; i < n; i++) {
      cerr << x[i] << " ";
    }
    cerr << endl;
    delete[] tmp1;
    delete[] tmp2;
    if (vector_norm_2(r1[1-j%2], n) < 1e-6) {
      cnt = j+1;
      goto out;
    }
  }
out:
  delete[] r1[0];
  delete[] r1[1];
  delete[] r2;
  delete[] p;
  return x;
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
  A = read_MatrixMarket("tests/sherman2.mtx", n);
  b = new double[n];
  for (int i = 0; i < n; i++) {
    b[i] = (1e-9*rand())/RAND_MAX;
  }
  x0 = new double[n];
  for (int i = 0; i < n; i++) {
    x0[i] = 0;
  }
  xm = bcgstab(A, b, x0, n, cnt);
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

