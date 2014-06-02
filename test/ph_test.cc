#include <phIO.h>
#include <apf.h>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv)
{
  int nodes, vars;
  ph_read_params(argv[1], argv[2], &nodes, &vars);
  printf("%d nodes, %d vars\n", nodes, vars);
  double* field;
  ph_read_field(argv[1], argv[2], &field);
  printf("first three values: %f %f %f\n", field[0], field[1], field[2]);
  ph_write_field("out", argv[2], field, nodes, vars, 0);
  free(field);
}
