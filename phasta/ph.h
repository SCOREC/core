#ifndef PH_H
#define PH_H

namespace ph {

void fail(const char* format, ...) __attribute__((noreturn,format(printf,1,2)));

}

#endif
