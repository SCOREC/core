#ifndef PH_IO_H
#define PH_IO_H

#ifdef __cplusplus
extern "C" {
#endif

enum chefio_file { CHEF_GEOMBC, CHEF_RESTART, NUM_CHEF_FILES};

int chefio_file_idx(const char* field);

#ifdef __INTEL_COMPILER
typedef size_t chefioTime;
#else
typedef struct timespec chefioTime;
#endif
void chefio_time(chefioTime* t);
size_t chefio_time_diff(chefioTime* start, chefioTime* end);
void chefio_setfile(int t);
void chefio_addOpenTime(size_t t);
void chefio_addCloseTime(size_t t);
void chefio_initStats();
void chefio_printStats();

void ph_write_preamble(FILE* f);
void ph_write_header(FILE* f, const char* name, size_t bytes,
    int nparam, int* params);

void ph_write_doubles(FILE* f, const char* name, double* data,
    size_t n, int nparam, int* params);
void ph_write_ints(FILE* f, const char* name, int* data,
    size_t n, int nparam, int* params);

/**
 * @brief determines if bytes read from the need to be 
 *        swapped to account for endianness
 * @return 1 if swapping is required, 0 otherwise
 */
int ph_should_swap(FILE* f);


/**
 *  @brief read a field
 *  @return 1 if no data block was read, 2 if data block read, 0 otherwise
 */
int ph_read_field(FILE* f, const char* field, int swap, 
    double** data, int* nodes, int* vars, int* step, char* hname);
void ph_write_field(FILE* f, const char* field, double* data,
    int nodes, int vars, int step);

#ifdef __cplusplus
}
#endif

#endif
