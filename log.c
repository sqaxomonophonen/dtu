#include "log.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

static void log_va(const char* fmt, va_list args) {
	FILE* out = stderr;
	vfprintf(out, fmt, args);
	fprintf(out, "\n");
	fflush(out);
}


void warnf(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_va(fmt, args);
	va_end(args);
}


void panicf(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_va(fmt, args);
	va_end(args);
	exit(EXIT_FAILURE);
}

