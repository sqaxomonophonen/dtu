#ifndef __A_H__
#define __A_H__

// a for assert, argh, abort, abandon ye all hope, and so on

void arghf(const char* fmt, ...) __attribute__((noreturn)) __attribute__((format (printf, 1, 2)));

#define ASSERT(cond) do { if(!(cond)) { arghf("ASSERT(%s) failed in %s() in %s:%d\n", #cond, __func__, __FILE__, __LINE__); } } while(0)

#define AN(expr) do { ASSERT((expr) != 0); } while(0)
#define AZ(expr) do { ASSERT((expr) == 0); } while(0)

#endif//__A_H__
