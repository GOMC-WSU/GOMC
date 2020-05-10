#ifndef ABORT_UNLESS_H
#define ABORT_UNLESS_H

#include <cstdio>
#include <cstdlib>

/*
 * Unlike the assert() macro, abort_unless() is not disabled by
 * defining the preprocessor macro NDEBUG.
 */

#define abort_unless(expr) do {                                         \
    if (!(expr)) {                                                      \
      std::fprintf(stderr, "%s:%u (%s): Assertion `%s' failed.\n",      \
                   __FILE__, __LINE__, __func__, __STRING(expr));       \
      std::fflush(stderr);                                              \
      abort();                                                          \
    }                                                                   \
  } while (0)

#endif
