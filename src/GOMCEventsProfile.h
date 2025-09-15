/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GOMCEVENTSPROFILE_H
#define GOMCEVENTSPROFILE_H

#include <map>

// Use the GOMC_PROFILE_EVENT macro within GOMCEventsProfileDef.h
// to define these event ID tags and string names.  The IDs generate
// a set of enumerated types and corresponding array of string names.
//
struct GomcProfileEvent {
  typedef enum {
#define GOMC_PROFILE_EVENT(a, b) a,
#include "GOMCEventsProfileDef.h"
#undef GOMC_PROFILE_EVENT
    EventsCount
  } Event;
};

char const *const GomcProfileEventStr[] = {
#define GOMC_PROFILE_EVENT(a, b) b,
#include "GOMCEventsProfileDef.h"
#undef GOMC_PROFILE_EVENT
    0};

#undef GOMC_PROFILE_START
#undef GOMC_PROFILE_STOP
#undef GOMC_REGISTER_EVENT
#undef GOMC_EVENT_START
#undef GOMC_EVENT_START_EX
#undef GOMC_EVENT_STOP
#undef GOMC_EVENT_RANGE
#undef GOMC_EVENT_RANGE_2

//
// Enable NVTX instrumentation for nvvp or Nsight profiling
// GOMC_CUDA build by defining GOMC_NVTX_ENABLED in CMake
//
#if defined(GOMC_CUDA) && defined(GOMC_NVTX_ENABLED)

#if CUDART_VERSION >= 10000
// #include
//</opt/nvidia/nsight-systems/2020.4.3/target-linux-x64/nvtx/include/nvtx3/nvToolsExt.h>
//// CUDA >= 10 has NVTX V3+
#include <nvtx3/nvToolsExt.h> // CUDA >= 10 has NVTX V3+
#else
#error NVTXv3 requires CUDA 10.0 or greater
// #include <nvToolsExt.h>        // CUDA < 10 has NVTX V2
#endif
#include <cuda_profiler_api.h>

// start profiling
#define GOMC_PROFILE_START()                                                   \
  do {                                                                         \
    cudaProfilerStart();                                                       \
  } while (0) // must terminate with semi-colon

// stop profiling
#define GOMC_PROFILE_STOP()                                                    \
  do {                                                                         \
    cudaProfilerStop();                                                        \
  } while (0) // must terminate with semi-colon

// C++ note: declaring const variables implies static (internal) linkage,
// and you have to explicitly specify "extern" to get external linkage.
const uint32_t GOMC_nvtx_colors[] = {
    0x0000ff00, 0x000000ff, 0x00ffff00, 0x00ff00ff, 0x0000ffff, 0x00ff0000,
    0x00006600, 0x00663300, 0x00000000, 0x007300e6, 0x00ff8c00,
};
const int GOMC_nvtx_colors_len = sizeof(GOMC_nvtx_colors) / sizeof(uint32_t);

#define GOMC_REGISTER_EVENT(name, cid)                                         \
  do {                                                                         \
  } while (0) // must terminate with semi-colon

// start recording an event
#define GOMC_EVENT_START(eon, id)                                              \
  do {                                                                         \
    if (eon) {                                                                 \
      int color_id = id;                                                       \
      color_id = color_id % GOMC_nvtx_colors_len;                              \
      nvtxEventAttributes_t eventAttrib = {0};                                 \
      eventAttrib.version = NVTX_VERSION;                                      \
      eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;                        \
      eventAttrib.colorType = NVTX_COLOR_ARGB;                                 \
      eventAttrib.color = GOMC_nvtx_colors[color_id];                          \
      eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;                       \
      eventAttrib.message.ascii = GomcProfileEventStr[id];                     \
      nvtxRangePushEx(&eventAttrib);                                           \
    }                                                                          \
  } while (0) // must terminate with semi-colon

// start recording an event
#define GOMC_EVENT_START_EX(eon, id, str)                                      \
  do {                                                                         \
    if (eon) {                                                                 \
      int color_id = id;                                                       \
      color_id = color_id % GOMC_nvtx_colors_len;                              \
      nvtxEventAttributes_t eventAttrib = {0};                                 \
      eventAttrib.version = NVTX_VERSION;                                      \
      eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;                        \
      eventAttrib.colorType = NVTX_COLOR_ARGB;                                 \
      eventAttrib.color = GOMC_nvtx_colors[color_id];                          \
      eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;                       \
      eventAttrib.message.ascii = str;                                         \
      nvtxRangePushEx(&eventAttrib);                                           \
    }                                                                          \
  } while (0) // must terminate with semi-colon

// stop recording an event
#define GOMC_EVENT_STOP(eon, id)                                               \
  do {                                                                         \
    if (eon) {                                                                 \
      nvtxRangePop();                                                          \
    }                                                                          \
  } while (0) // must terminate with semi-colon

// embed event recording in class to automatically pop when destroyed
class GOMC_NVTX_Tracer {
protected:
  int evon; // is event on?
public:
  GOMC_NVTX_Tracer(int eon, int id = 0) : evon(eon) {
    GOMC_EVENT_START(eon, id);
  }
  ~GOMC_NVTX_Tracer() { GOMC_EVENT_STOP(evon, 0); }
};

// call GOMC_EVENT_RANGE at beginning of function to push event recording
// destructor is automatically called on return to pop event recording
#define GOMC_EVENT_RANGE(eon, id) GOMC_NVTX_Tracer namd_nvtx_tracer(eon, id)
// must terminate with semi-colon

#if defined(GOMC_PROFILE_EVENT_LAYER_2)
#define GOMC_EVENT_RANGE_2(eon, id) GOMC_EVENT_RANGE(eon, id)
#else
#define GOMC_EVENT_RANGE_2(eon, id)                                            \
  do {                                                                         \
  } while (0) // must terminate with semi-colon
#endif

#else
//
// Otherwise all profiling macros become no-ops.
//
#define GOMC_PROFILE_START()                                                   \
  do {                                                                         \
  } while (0) // must terminate with semi-colon

#define GOMC_PROFILE_STOP()                                                    \
  do {                                                                         \
  } while (0) // must terminate with semi-colon

#define GOMC_REGISTER_EVENT(name, cid)                                         \
  do {                                                                         \
  } while (0) // must terminate with semi-colon

#define GOMC_EVENT_START(eon, id)                                              \
  do {                                                                         \
  } while (0) // must terminate with semi-colon

#define GOMC_EVENT_START_EX(eon, id, str)                                      \
  do {                                                                         \
  } while (0) // must terminate with semi-colon

#define GOMC_EVENT_STOP(eon, id)                                               \
  do {                                                                         \
  } while (0) // must terminate with semi-colon

#define GOMC_EVENT_RANGE(eon, id)                                              \
  do {                                                                         \
  } while (0) // must terminate with semi-colon

#define GOMC_EVENT_RANGE_2(eon, id)                                            \
  do {                                                                         \
  } while (0) // must terminate with semi-colon

#endif // GOMC_CUDA && GOMC_NVTX_ENABLED

#endif /* GOMCEVENTSPROFILE_H */
