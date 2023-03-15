#include "FDM.h"

double Functions::func(double x, double y) {
   return -2 + x * x;
}

double Functions::u_g(double x, double y) {
   return x * x;
}

double Functions::u_betta(double x, double y) {
   return x * x;
}

double Functions::u_real(double x, double y) {
   return x * x;
}

/* Конструктор */
FDM::FDM() {
   hx = 0; hy = 0;
   nx = 0; ny = 0;
   N = 0;
   kx = 0; ky = 0; 
   hx_i = 0; hy_j = 0;
   hx0 = 0; hy0 = 0;
   nodes_x = vector<double>();
   nodes_y = vector<double>();
   nodes_xT = vector<double>();
   nodes_yT = vector<double>();
   u_real = vector<double>();
   u_answer = vector<double>();
   F = vector<double>();
   A = vector<vector<double>>();
   first_bc = vector<int>();
   thrid_bc = vector<int>();
   gamma = 0; lambda = 0; beta = 0;
}

/* Деструктор */
FDM::~FDM() {
   nodes_x.clear();
   nodes_y.clear();
   nodes_xT.clear();
   nodes_yT.clear();
   A.clear();
   F.clear();
   u_real.clear();
   u_answer.clear();
}

/* РАВНОМЕРНАЯ СЕТКА
 * На ввод подаются коэффициенты гамма и лямбда,
 * количество узлов сетки по х и по у,
 * х и у координаты узлов, в которых
 * имеет границы расчетная сетка, 
 * имеющая форму буквы Т.
 * Рассчитывается сетка с фиктивными узлами -
 * то есть сетка прямоугольной формы, для построения матрицы.
 * Рассчитывается вектор правой части.
 */
void FDM::read_data() {
   Functions B;
   ifstream fin_data("data.txt");

   fin_data >> lambda >> gamma;
   fin_data >> nx >> ny;

   nodes_xT.resize(4);
   nodes_yT.resize(3);
   nodes_x.resize(nx);
   nodes_y.resize(ny);
   
   N = nx * ny;

   u_answer.resize(N);
   u_real.resize(N);
   F.resize(N);

   for (int i = 0; i < 4; i++) {
      fin_data >> nodes_xT[i];
   }
   for (int i = 0; i < 3; i++) {
      fin_data >> nodes_yT[i];
   }

   fin_data.close();

   sort(nodes_xT.begin(), nodes_xT.end());
   sort(nodes_yT.begin(), nodes_yT.end());

   hx = (nodes_xT[3] - nodes_xT[0]) / (nx - 1.0);
   hy = (nodes_yT[2] - nodes_yT[0]) / (ny - 1.0);

   for (int i = 0; i < nx; i++) nodes_x[i] = nodes_xT[0] + i * hx;
   for (int i = 0; i < ny; i++) nodes_y[i] = nodes_yT[0] + i * hy;

   for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
         F[i + nx * j] = B.func(nodes_x[i], nodes_y[j]);
      }
   }
}

/* РАВНОМЕРНАЯ СЕТКА
 * Заполняем матрицу значениями (вместе с
 * фиктивными узлами), затем производим учет
 * фиктивных узлов - обнуляем соответствующие элементы
 * вектора u.
 */
void FDM::build_matrix() {
   A.resize(5);
   for (int i = 0; i < 5; i++) A[i].resize(N);

   double Lx = -lambda / (hx * hx);
   double Ly = -lambda / (hy * hy);
   double D = 2.0 * lambda * (1.0 / (hx * hx) + 1.0 / (hy * hy)) + gamma;

   for (int i = nx; i < N; i++) A[0][i] = Ly;
   for (int i = 1; i < N; i++) A[1][i] = Lx;
   for (int i = 0; i < N; i++) A[2][i] = D;
   for (int i = 0; i < N - 1; i++) A[3][i] = Lx;
   for (int i = 0; i < N - nx; i++) A[4][i] = Ly;

   for (int i = 0; nodes_y[i] < nodes_yT[1]; i++) {
      for (int j = 0; j < nx; j++) {
         if (nodes_x[j] < nodes_xT[1] || nodes_x[j] > nodes_xT[2]) {
            A[0][j + i * nx] = A[1][j + i * nx] = A[3][j + i * nx] = A[4][j + i * nx] = 0.0;
            A[2][j + i * nx] = 1.0;
            F[j + i * nx] = 0.0;
         }
      }
   } 
}

/* НЕРАВНОМЕРНАЯ СЕТКА
 * На ввод подаются коэффициенты гамма и лямбда,
 * первый шаг по x и по y, множители для 
 * последующих шагов, количество узлов сетки по х и по у,
 * х и у координаты узлов, в которых
 * имеет границы расчетная сетка, имеющая форму буквы Т.
 * Рассчитывается сетка с фиктивными узлами -
 * то есть сетка прямоугольной формы, для построения матрицы.
 * Рассчитывается вектор правой части.
 */
void FDM::uneven_grid_read_data() {
   Functions B;
   ifstream fin_un_data("data_uneven.txt");

   fin_un_data >> lambda >> gamma;
   fin_un_data >> hx0 >> kx >> nx;
   fin_un_data >> hy0 >> ky >> ny;
   
   nodes_xT.resize(4);
   nodes_yT.resize(3);
   nodes_x.resize(nx);
   nodes_y.resize(ny);

   N = nx * ny;

   u_answer.resize(N);
   u_real.resize(N);
   F.resize(N);

   for (int i = 0; i < 4; i++) {
      fin_un_data >> nodes_xT[i];
   }
   for (int i = 0; i < 3; i++) {
      fin_un_data >> nodes_yT[i];
   }

   fin_un_data.close();

   nodes_x[0] = nodes_xT[0];
   for (int i = 1; i < nx; i++) {
      hx_i = hx0 * pow(kx, i - 1);
      nodes_x[i] = nodes_x[i - 1] + hx_i;
   }
   nodes_y[0] = nodes_yT[0];
   for (int i = 1; i < ny; i++) {
      hy_j = hy0 * pow(ky, i - 1);
      nodes_y[i] = nodes_y[i - 1] + hy_j;
   }

   for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
         F[i + nx * j] = B.func(nodes_x[i], nodes_y[j]);
      }
   }
}

/* НЕРАВНОМЕРНАЯ СЕТКА
 * Заполняем матрицу значениями (вместе с
 * фиктивными узлами), затем производим учет
 * фиктивных узлов - обнуляем соответствующие элементы
 * вектора u.
 */
void FDM::build_uneven_matrix() {
   double tmp = 0.0;
   hx_i = hx0; hy_j = hy0;

   A.resize(5);
   for (int i = 0; i < 5; i++) A[i].resize(N);

   for (int i = nx; i < N; i++) {
      tmp = hy_j; // hy_(j - 1)
      hy_j = hy0 * pow(ky, i / nx); // hy_j
      A[0][i] = -2 * lambda * (1.0 / (tmp * (hy_j + tmp)));
   }

   hx_i = hx0; hy_j = hy0;
   for (int i = 1; i < N; i++) {
      tmp = hx_i;  // hx_(i - 1)
      hx_i = hx0 * pow(kx, i % nx);  // hx_i 
      A[1][i] = -2 * lambda * (1.0 / (tmp * (hx_i + tmp)));
   }

   hx_i = hx0; hy_j = hy0;
   for (int i = 0; i < N; i++) {
      tmp = hx_i;
      hx_i = hx0 * pow(kx, i % nx);  // hx_i  
      A[2][i] = 2.0 * lambda * (1.0 / (hx_i * tmp)) + gamma;
      tmp = hy_j;
      hy_j = hy0 * pow(ky, i / nx); // hy_j
      A[2][i] += 2.0 * lambda * (1.0 / (hy_j * tmp));
   }

   hx_i = hx0; hy_j = hy0;
   for (int i = 0; i < N - 1; i++) {
      tmp = hx_i;  // hx_(i - 1)
      hx_i = hx0 * pow(kx, i % nx);  // hx_i 
      A[3][i] = -2 * lambda * (1.0 / (hx_i * (hx_i + tmp)));
   }

   hx_i = hx0; hy_j = hy0;
   for (int i = 0; i < N - nx; i++) {
      tmp = hy_j;
      hy_j = hy0 * pow(ky, i / nx);
      A[4][i] = -2 * lambda * (1.0 / (hy_j * (hy_j + tmp)));
   }

   for (int i = 0; nodes_y[i] < nodes_yT[1]; i++) {
      for (int j = 0; j < nx; j++) {
         if (nodes_x[j] < nodes_xT[1] || nodes_x[j] > nodes_xT[2]) {
            A[0][j + i * nx] = A[1][j + i * nx] = A[3][j + i * nx] = A[4][j + i * nx] = 0.0;
            A[2][j + i * nx] = 1.0;
            F[j + i * nx] = 0.0;
         }
      }
   }
}


/*
 * Чтение краевых условий из файла,
 * представлены краевые условия первого и третьего рода.
 */
void FDM::read_BC() {
   ifstream fin_boundary("boundary.txt");
   int count = 0;  // колличество граничных узлов для КУ
   int type = 0;
   
   fin_boundary >> beta;

   fin_boundary >> type >> count;
   if (type == 1) {
      first_bc.resize(count);
      for (int i = 0; i < count; i++) {
         fin_boundary >> first_bc[i];
      }
      sort(first_bc.begin(), first_bc.end());
      type = 0;
      count = 0;
   }

   fin_boundary >> type >> count;
   if (type == 3) {
      thrid_bc.resize(count);
      for (int i = 0; i < count; i++) {
         fin_boundary >> thrid_bc[i];
      }
      sort(thrid_bc.begin(), thrid_bc.end());
   }

   fin_boundary.close();
}

/*
 * Производится учет краевых условий 
 */
void FDM::boundary_condition() {
   double Hb = 0.0, H = 0.0, x = 0.0, y = 0.0;
   Functions func;

   // Третье краевое
   for (int i = 0; i < thrid_bc.size(); i++) {
      x = nodes_x[thrid_bc[i] % nx];
      y = nodes_y[thrid_bc[i] / nx];
      switch (direction_n(x, y)) {
         case 1: {
            Hb = beta + (lambda / hx);
            H = - (lambda / hx);
            A[0][thrid_bc[i]] = A[1][thrid_bc[i]] = A[4][thrid_bc[i]] = 0.0;
            A[2][thrid_bc[i]] = Hb;
            A[3][thrid_bc[i]] = H;
            break;
         }
         case 2: {
            Hb = beta + (lambda / hx);
            H = -(lambda / hx);
            A[0][thrid_bc[i]] = A[3][thrid_bc[i]] = A[4][thrid_bc[i]] = 0.0;
            A[2][thrid_bc[i]] = Hb;
            A[1][thrid_bc[i]] = H;
            break;
         }
         case 3: {
            Hb = beta + (lambda / hy);
            H = -(lambda / hy);
            A[0][thrid_bc[i]] = A[1][thrid_bc[i]] = A[3][thrid_bc[i]] = 0.0;
            A[2][thrid_bc[i]] = Hb;
            A[4][thrid_bc[i]] = H;
            break;
         }
         case 4: {
            Hb = beta + (lambda / hy);
            H = -(lambda / hy);
            A[3][thrid_bc[i]] = A[1][thrid_bc[i]] = A[4][thrid_bc[i]] = 0.0;
            A[2][thrid_bc[i]] = Hb;
            A[0][thrid_bc[i]] = H;
            break;
         }
         default:
            break;
      }
      F[thrid_bc[i]] = beta * func.u_betta(x, y);
   }

   // Первое краевое
   for (int i = 0; i < first_bc.size(); i++) {
      x = nodes_x[first_bc[i] % nx];
      y = nodes_y[first_bc[i] / nx];
      A[0][first_bc[i]] = A[1][first_bc[i]] = A[3][first_bc[i]] = A[4][first_bc[i]] = 0.0;
      A[2][first_bc[i]] = 1.0;
      F[first_bc[i]] = func.u_g(x, y);
   }
}

/*
 * Проверка направления ветора нормали
 */
int FDM::direction_n(double x, double y) {
   int dir = -1;
   if ((x == nodes_xT[0] || x == nodes_xT[1]) && y != nodes_yT[2])
      dir = 1; // влево
   else if ((x == nodes_xT[2] || x == nodes_xT[3]) && y != nodes_yT[2])
      dir = 2; // вправо
   else if (y == nodes_yT[0] || y == nodes_yT[1])
      dir = 3; // вниз
   else if (y == nodes_yT[2])
      dir = 4; // вверх
   return dir;
}

/* 
 * Вывод результатов производится в файл
 * типа .csv в виде таблицы, в этот же файл выводятся
 * значения искомой функции во всех точках
 * и абсолютная погрешность.
 */
void FDM::print_result() {
   double error = 0.0;
   Functions B;

   ofstream fout_answ("out_answ.csv");
   fout_answ.precision(16);

   fout_answ << "answer" << ";" << "real" << ";" << "error" << ";" << endl;

   for (int i = 0; i < ny; i++) {
      for (int j = 0; j < nx; j++) {
         fout_answ << u_answer[j + i * nx] << ";";
         if (nodes_x[j] >= nodes_xT[1] && nodes_x[j] <= nodes_xT[2] || 
             nodes_y[i] >= nodes_yT[1]) {
            u_real[j + i * nx] = B.u_real(nodes_x[j], nodes_y[i]);
         }
         fout_answ << u_real[j + i * nx] << ";";
         fout_answ << u_real[j + i * nx] - u_answer[j + i * nx] << ";" << endl;
      }
   }
   fout_answ.close();
}

