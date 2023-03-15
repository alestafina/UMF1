#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

class FDM {
protected: 
   double hx, hy;                       // Шаги по x и y
   double hx0, hy0;                     // Начальный шаг по x и по y для неравномерной сетки
   double kx, ky, hx_i, hy_j;           // Множитель шага, шаг по x и шаг по y 
   int nx, ny, N;                       // Количество узлов и размерность матрицы
   vector<double> nodes_x, nodes_y;     // Координаты узлов по x, по y
   vector<double> nodes_xT, nodes_yT;   // Координаты узлов сетки формы Т
   vector<double> u_real, u_answer;     // Настоящее значение u и полученное значение u
   vector<vector<double>> A;            // Матрица А в диагональном формате
   vector<double> F;                    // Вектор правой части
   vector<int> first_bc, thrid_bc;      // Номера узлов с соответствующими КУ
   double gamma, lambda, beta;          // Коэффициенты гамма и лямбда
public:                                 
   FDM();                               // Конструктор по умолчанию
   ~FDM();                              // Деструктор 
                                        
   void read_data();                    // Считывание данных
   void build_matrix();                 // Генерация матрицы
   void read_BC();                      // Считывание краевых условий

   void uneven_grid_read_data();        // Чтение данных нерегулярной сетки
   void build_uneven_matrix();          // Генерация матрицы для неравномерной сетки

   void boundary_condition();           // Учет первых и третьих краевых условий
   int direction_n(double x, double y); // Направление нормали (грань)
   void print_result();
};

class Functions {
public:
   double func(double x, double y);
   double u_g(double x, double y);
   double u_betta(double x, double y);
   double u_real(double x, double y);
};

class Solve : public FDM {
private: 
   double eps, w, residual, norm;
   int max_iter;
public:
   Solve();
   ~Solve();

   void iteration_step(vector<double> &x0, vector<double> &x1, double w);
   void seidel();
};
