#ifndef __SCRATCH_H__
#define __SCRATCH_H__

#include <strings.h>

struct scratch {
	size_t capacity;
	size_t top;
	void* ptr;
};

struct scratch_lap {
	struct scratch* scratch;
	size_t top;
};


#define SCRATCH_MIN_CAPACITY (1<<20)

static void _scratch_ensure(struct scratch* scratch)
{
	if(scratch->top > scratch->capacity) {
		if(scratch->capacity < SCRATCH_MIN_CAPACITY) {
			scratch->capacity = scratch->top > SCRATCH_MIN_CAPACITY ? scratch->top : SCRATCH_MIN_CAPACITY;
		} else {
			scratch->capacity <<= 1;
		}
		//printf("DEBUG -- new capacity: %zd\n", scratch->capacity);
		scratch->ptr = realloc(scratch->ptr, scratch->capacity);
	}
}

static inline void scratch_init(struct scratch* scratch)
{
	bzero(scratch, sizeof(struct scratch));
}

static inline size_t scratch_align(struct scratch* scratch, size_t element_sz)
{
	size_t rem = scratch->top % element_sz;
	if(rem > 0) scratch->top += (element_sz - rem);
	return scratch->top;
}

static inline void* scratch_get(struct scratch* scratch, size_t offset)
{
	return (void*)((char*)scratch->ptr + offset);
}

static inline void* scratch_alloc(struct scratch* scratch, size_t element_sz)
{
	scratch->top += element_sz;
	_scratch_ensure(scratch);
	return scratch_get(scratch, scratch->top - element_sz);
}


static inline void scratch_reset(struct scratch* scratch)
{
	scratch->top = 0;
}

#endif//__SCRATCH_H__
