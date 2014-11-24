double** ilu(double** A, int n) {
  double** B = new double*[n];
  for (int i = 0; i < n; i++) {
    B[i] = new double[n];
    for (int j = 0; j < n; j++) {
      B[i][j] = A[i][j];
    }
  }
  for (int i = 1; i < n; i++) {
    for (int k = 0; k < i; k++) {
      if (fabs(A[i][k]) < 1e-15) continue;
      if (fabs(B[k][k]) < 1e-15) {
        for (int j = 0; j < n; j++) {
          delete[] B[j];
        }
        delete[] B;
        return NULL;
      }
      B[i][k] = B[i][k]/B[k][k];
      for (int j = k+1; j < n; j++) {
        if (fabs(A[i][j]) < 1e-15) continue;
        B[i][j] = B[i][j] - B[i][k]*B[k][j];
      }
    }
  }
  return B;
}
