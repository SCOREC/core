#ifndef SCOREC_PCU_PCU_DEFINES_H
#define SCOREC_PCU_PCU_DEFINES_H

#define PCU_SUCCESS 0
#define PCU_FAILURE -1
#ifdef __GNUC__
#define PCU_FORMAT_ATTRIBUTE(...) \
  __attribute__((format(printf, ##__VA_ARGS__)));
#else
#define PCU_FORMAT_ATTRIBUTE(format, ...)
#endif

#ifdef __cplusplus
extern "C"{
#endif

struct PCU_t {
    void* ptr;
};

#ifdef __cplusplus
}
#endif

#endif // SCOREC_PCU_PCU_DEFINES_H