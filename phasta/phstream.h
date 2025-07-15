#ifndef PHSTREAM_H_
#define PHSTREAM_H_
#include <stdio.h>
#include <PCU.h>


/** \file phstream.h
    \brief The data stream API
    \remark This api supports the creation and destruction and
            manipulation of data streams as used in chefStream.cc and
            <a href=https://github.com/PHASTA/phastaChef>phastaChef</a>.
*/

typedef struct RStream* rstream;
typedef struct GRStream* grstream;
/** @brief make restart stream */
rstream makeRStream(pcu::PCU *PCUObj);
/** @brief clear restart stream */
void clearRStream(rstream rs, pcu::PCU *PCUObj);
/** @brief detach output stream */
void destroyRStream(rstream rs, pcu::PCU *PCUObj);

/** @brief make geom-restart stream */
grstream makeGRStream(pcu::PCU *PCUObj);
/** @brief clear geom-restart stream */
void clearGRStream(grstream grs, pcu::PCU *PCUObj);
/** @brief destroy geom-restart stream */
void destroyGRStream(grstream grs, pcu::PCU *PCUObj);

/** @brief open restart stream for reading*/
FILE* openRStreamRead(rstream rs, pcu::PCU *PCUObj);
/** @brief open restart stream for writing*/
FILE* openRStreamWrite(rstream rs, pcu::PCU *PCUObj);

/** @brief open named stream in geom-restart stream for reading*/
FILE* openGRStreamRead(grstream grs, const char* named, pcu::PCU *PCUObj);
/** @brief open named stream in geom-restart stream for writing*/
FILE* openGRStreamWrite(grstream grs, const char* named, pcu::PCU *PCUObj);

/** @brief dev function */
void attachRStream(grstream grs, rstream rs, pcu::PCU *PCUObj);
#endif 
