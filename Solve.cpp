#include "FDM.h"

/*
 * Решение СЛАУ методом Гаусса-Зейделя
 */
Solve::Solve() {
   eps = 1e-16;
   w = 1.66;
   residual = 1.0;
   norm = 0.0;
   max_iter = 1000;
   u_answer = vector<double>();
}

Solve::~Solve() {
   u_answer.clear();
}

void Solve::iteration_step(vector<double> &x0, vector<double> &x1, double w) {
   int ind[5] = { -nx, -1, 0, 1, nx };
   double sum = 0.0;
   int i, j, k;

   for (i = 0; i < N; i++) {
      sum = 0.0, residual = 0.0;
      for (j = 0; j < 5; j++) {
         k = ind[j] + i;
         if (k >= 0 && k < N)
            sum += A[j][i] * x0[k];
      }
      x1[i] = x0[i] + w * (F[i] - sum) / A[2][i];
      residual += (sum - F[i]) * (sum - F[i]);
   }
   residual = sqrt(residual / norm);
}

void Solve::seidel() {
   u_answer.resize(N);
   for (int i = 0; i < N; i++) norm += F[i] * F[i];
   for (int k = 0; k < max_iter && residual > eps; k++)
      iteration_step(u_answer, u_answer, w);
}
