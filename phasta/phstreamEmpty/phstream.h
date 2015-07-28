#ifndef PHSTREAM_EMPTY_H_
#define PHSTREAM_EMPTY_H_
struct RStream;
struct GRStream;
/** @brief open restart stream for reading*/
FILE* openRStreamRead(RStream* rs);
/** @brief open named stream in geom-restart stream for writing*/
FILE* openGRStreamWrite(GRStream* grs, const char* named);
#endif 
