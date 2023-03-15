#include "FDM.h"

int main() {
   Solve matrix;
   matrix.read_data();
   matrix.build_matrix();
   matrix.read_BC();
   matrix.boundary_condition();
   matrix.seidel();
   matrix.print_result();
   return 0;
}
