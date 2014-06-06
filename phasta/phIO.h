#ifndef PH_IO_H
#define PH_IO_H

#ifdef __cplusplus
extern "C" {
#endif

void ph_write_preamble(FILE* f);
void ph_write_header(FILE* f, const char* name, size_t bytes,
    int nparam, int* params);

void ph_write_doubles(FILE* f, const char* name, double* data,
    size_t n, int nparam, int* params);
void ph_write_ints(FILE* f, const char* name, int* data,
    size_t n, int nparam, int* params);

void ph_read_field(const char* file, const char* field, double** data,
    int* nodes, int* vars, int* step);
void ph_write_field(FILE* f, const char* field, double* data,
    int nodes, int vars, int step);

#ifdef __cplusplus
}
#endif

#endif
