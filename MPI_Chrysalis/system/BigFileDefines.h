#ifndef BIG_FILE_DEFINES_H
#define BIG_FILE_DEFINES_H

#ifdef __linux
     #ifndef _LARGEFILE64_SOURCE
     #define _LARGEFILE64_SOURCE
     #endif
     #define _FILE_OFFSET_BITS 64
     #define STAT_NAME stat64
#else
     #define STAT_NAME stat
#endif

#endif //BIG_FILE_DEFINES_H
