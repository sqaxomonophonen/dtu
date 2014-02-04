#ifndef __PSYS_H__
#define __PSYS_H__

#include <stdint.h>

#include "scratch.h"

struct aabb {
	float x0, y0;
	float x1, y1;
};

struct solid_op {
	enum {
		SOP_FRACTURE
	} type;
	union {
		struct {
			int foo, bar;
		} fracture;
	};
	struct solid_op* next;
};

struct solid {
	float px, py; // position
	float vx, vy; // velocity
	float fx, fy; // forces
	float r; // rotation
	float vr; // rotational velocity
	float fvr; // torque

	float m; // mass
	// center of mass (x,y) = (zmx/m, zmy/m)? (to make changes easy)
	float zmx, zmy; // (zmx/m, zmy/m) is center of mass
	float I; // moment of inertia

	struct aabb aabb;

	// bitmap dimensions
	int b_width, b_height;

	// bitmap color
	uint32_t* b_rgba;
	GLuint gl_rgba;

	// bitmap cell types, 0=unoccupied, 1=particle collision halo, 2=solid type #1, etc?
	uint8_t* b_type;

	// scratch buffer? make b_type wider (32bit)? could possibly encode a
	// lot of stuff. though, it probably isn't necessary to pack stuff in a
	// wide word (if L1/L2 caches are my concern; the CPU can have multiple
	// memory ranges in cache at the same time)

	// linked list of running operations, e.g. ongoing fractures
	struct solid_op* ops;

	uint32_t dirty_flags;


	// derived values
	float inv_m; // inverse mass
	float inv_I; // inverse moment of inertia;
	float cx, cy; // center of mass
	float tx_x0, tx_y0, tx_u, tx_v; // transform


	// next solid in linked list?
	struct solid* next;
};


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
	int occupied_buckets_count;

	int collisions;

	struct solid* solids;

	struct scratch scratch;
};

void psys_init(struct psys* ps);
void psys_step(struct psys* ps);
void psys_draw(struct psys* ps);


int solid_normal_at_world_point(struct solid* solid, float px, float py, float* nx, float* ny);

#endif//__PSYS_H__
