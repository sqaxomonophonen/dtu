#ifndef __LOG_H__
#define __LOG_H__

#define PRINTF_LIKE_FMT __attribute__((format(printf, 1, 2)))
#define NO_RETURN __attribute__((noreturn))

void warnf(const char* fmt, ...) PRINTF_LIKE_FMT;
void panicf(const char* fmt, ...) PRINTF_LIKE_FMT NO_RETURN;

#endif//__LOG_H__
