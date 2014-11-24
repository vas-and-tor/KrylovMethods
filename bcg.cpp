#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "matrix.h"
#include <cmath>
#include <cstdio>

using namespace std;

double * bcg(double** A, double* b, double* x0, int n, int & cnt) {
  double* r1[2] = { NULL         , new double[n] };
  double* r2[2] = { new double[n], new double[n] };
  double* p1[2] = { new double[n], new double[n] };
  double* p2[2] = { new double[n], new double[n] };
  double* x [2] = { new double[n], new double[n] };
  double** AT = new double*[n];
  bool converges = false;
  r1[0] = mult_matrix_to_vector(A, x0, n);
  for (int i = 0; i < n; i++) {
    p1[0][i] = p2[0][i] = r2[0][i] = r1[0][i] = b[i] - r1[0][i];
    AT[i] = new double[n];
    x[0][i] = x0[i];
  }
  if (mult_vector_to_vector(r1[0], r2[0], n) < 1e-9) {
    converges = true;
    goto out;
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      AT[i][j] = A[j][i];
    }
  }
  for (int j = 0; j < 300; j++) {
    double* tmp1 = mult_matrix_to_vector(A , p1[j%2], n);
    double* tmp2 = mult_matrix_to_vector(AT, p2[j%2], n);
    double alpha = mult_vector_to_vector(r1[j%2], r2[j%2], n)/mult_vector_to_vector(tmp1, p2[j%2], n);
    double betha;
    for (int i = 0; i < n; i++) {
      r1[1-j%2][i] = r1[j%2][i];
      r2[1-j%2][i] = r2[j%2][i];
    }
    tmp1 = mult_vector_to_scalar_inplace(tmp1, -alpha, n);
    tmp2 = mult_vector_to_scalar_inplace(tmp2, -alpha, n);
    r1[1-j%2] = add_vector_to_vector_inplace(r1[1-j%2], tmp1, n);
    r2[1-j%2] = add_vector_to_vector_inplace(r2[1-j%2], tmp2, n);
    betha = mult_vector_to_vector(r1[1-j%2], r2[1-j%2], n)/mult_vector_to_vector(r1[j%2], r2[j%2], n);
    for (int i = 0; i < n; i++) {
      x [1-j%2][i] = x [  j%2][i] + alpha*p1[j%2][i];
      p1[1-j%2][i] = r1[1-j%2][i] + betha*p1[j%2][i];
      p2[1-j%2][i] = r2[1-j%2][i] + betha*p2[j%2][i];
    }
    cerr << j+1 << " Answer: ";
    for (int i = 0; i < n; i++) {
      cerr << x[1-j%2][i] << " ";
    }
    cerr << endl;
    cerr << vector_norm_2(r1[1-j%2], n) << " " << fabs(betha) << endl;
    delete[] tmp1;
    delete[] tmp2;
    if (vector_norm_2(r1[1-j%2], n) < 1e-6/* || fabs(betha) < 1e-15*/) {
      converges = (vector_norm_2(r1[1-j%2], n) < 1e-6);
      cnt = j+1;
      goto out;
    }
  }
out:
  for (int i = 0; i < n; i++) {
    delete[] AT[i];
  }
  delete[] AT;
  for (int i = 0; i < 2; i++) {
    delete[] r1[i];
    delete[] r2[i];
    delete[] p1[i];
    delete[] p2[i];
  }
  delete x[1-cnt%2];
  double* xm = x[cnt%2];
  cnt = (converges ? cnt : -1);
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
    b[i] = (1e3*rand())/RAND_MAX;
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

