#ifndef PH_IO_H
#define PH_IO_H

#ifdef __cplusplus
extern "C" {
#endif

void ph_read_params(
    const char* file,
    const char* field,
    int* nodes, int* vars);
void ph_read_field(
    const char* file,
    const char* field,
    double** data);
void ph_write_field(
    const char* file,
    const char* field,
    double* data,
    int nodes,
    int vars,
    int step);

#ifdef __cplusplus
}
#endif

#endif
