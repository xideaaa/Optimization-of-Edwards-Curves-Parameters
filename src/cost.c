#include "cost.h"

#ifdef FF_COUNT_OPS
_Thread_local uint64_t FF_COUNT_MULS = 0;
_Thread_local uint64_t FF_COUNT_SQRS = 0;
_Thread_local uint64_t FF_COUNT_ADDS = 0;
#endif
