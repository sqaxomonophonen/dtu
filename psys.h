#ifndef __PSYS_H__
#define __PSYS_H__

#include <stdint.h>

struct ppair {
	int i, j;
};

struct particle {
	float px;
	float py;
	uint32_t cell_key; // doubles as occupancy marker?
	float vx;
	float vy;
	float vxe;
	float vye;
	float density;
};

struct particle_hash_bucket {
	struct particle particle;
	int insert_offset;
};

struct particle_hash {
	struct particle_hash_bucket* buckets;
};

struct psys {
	int particle_count;
	struct particle_hash hashes[2];
	int current_hash_index;

	int bucket_count_shift;

	int* occupied_buckets;
	int occupied_buckets_max;
	int occupied_buckets_count;

	int collisions;

	int ppair_count;
	struct ppair* ppairs;
};

void psys_init(struct psys* ps);
void psys_step(struct psys* ps);
void psys_draw(struct psys* ps);

#endif//__PSYS_H__
