#ifndef LION_PRINT_H
#define LION_PRINT_H
#include <stdio.h>
#include <stdarg.h>
#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief set the verbosity of output
 * \param lvl (in) 0 = off, increasing values (> 0) increase the amount of output
 */
void lion_set_verbosity(int lvl);

/**
 * \brief get the verbosity level for output
 */
int lion_get_verbosity();

/**
 * \brief set the stdout file stream destination
 * \param out (in) the file stream written to by lion_oprint
 */
void lion_set_stdout(FILE* out);

/**
 * \brief set the stderr file stream destination
 * \param err (in) the file stream written to by lion_eprint
 */
void lion_set_stderr(FILE* err);

/**
 * \brief fprintf(stdout,...) wrapper
 * \param lvl (in) the verbosity level of this message
 * \remark see printf for arguments
 */
int lion_oprint(int lvl, char const*, ...);

/**
 * \brief fprintf(stderr,...) wrapper
 * \param lvl (in) the verbosity level of this message
 * \remark see printf for arguments
 */
int lion_eprint(int lvl, char const*, ...);

/**
 * \brief variadic fprintf(stdout,...) wrapper
 * \param lvl (in) the verbosity level of this message
 * \remark see printf for arguments, used for nesting within variadic functions
 */
int lion_voprint(int lvl, char const*, va_list);

/**
 * \brief variadic fprintf(stderr,...) wrapper
 * \param lvl (in) the verbosity level of this message
 * \remark see printf for arguments, used for nesting within variadic functions
 */
int lion_veprint(int lvl, char const*, va_list);

#ifdef __cplusplus
}
#endif
#endif // header guard
