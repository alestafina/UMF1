#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

class FDM {
protected: 
   double hx, hy;                       // ���� �� x � y
   double hx0, hy0;                     // ��������� ��� �� x � �� y ��� ������������� �����
   double kx, ky, hx_i, hy_j;           // ��������� ����, ��� �� x � ��� �� y 
   int nx, ny, N;                       // ���������� ����� � ����������� �������
   vector<double> nodes_x, nodes_y;     // ���������� ����� �� x, �� y
   vector<double> nodes_xT, nodes_yT;   // ���������� ����� ����� ����� �
   vector<double> u_real, u_answer;     // ��������� �������� u � ���������� �������� u
   vector<vector<double>> A;            // ������� � � ������������ �������
   vector<double> F;                    // ������ ������ �����
   vector<int> first_bc, thrid_bc;      // ������ ����� � ���������������� ��
   double gamma, lambda, beta;          // ������������ ����� � ������
public:                                 
   FDM();                               // ����������� �� ���������
   ~FDM();                              // ���������� 
                                        
   void read_data();                    // ���������� ������
   void build_matrix();                 // ��������� �������
   void read_BC();                      // ���������� ������� �������

   void uneven_grid_read_data();        // ������ ������ ������������ �����
   void build_uneven_matrix();          // ��������� ������� ��� ������������� �����

   void boundary_condition();           // ���� ������ � ������� ������� �������
   int direction_n(double x, double y); // ����������� ������� (�����)
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
