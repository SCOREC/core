#ifndef PHSTREAM_H_
#define PHSTREAM_H_
#include <stdio.h>
typedef struct RStream* rstream;
typedef struct GRStream* grstream;
/** @brief make restart stream */
rstream makeRStream();
/** @brief clear restart stream */
void clearRStream(rstream rs);
/** @brief detach output stream */
void destroyRStream(rstream rs);

/** @brief make geom-restart stream */
grstream makeGRStream();
/** @brief clear geom-restart stream */
void clearGRStream(grstream grs);
/** @brief destroy geom-restart stream */
void destroyGRStream(grstream grs);

/** @brief open restart stream for reading*/
FILE* openRStreamRead(rstream rs);
/** @brief open restart stream for writing*/
FILE* openRStreamWrite(rstream rs);

/** @brief open named stream in geom-restart stream for reading*/
FILE* openGRStreamRead(grstream grs, const char* named);
/** @brief open named stream in geom-restart stream for writing*/
FILE* openGRStreamWrite(grstream grs, const char* named);

/** @brief dev function */
void attachRStream(grstream grs, rstream rs);
#endif 
