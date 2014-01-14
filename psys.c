#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>

#include <glew.h>

#include "psys.h"

// XXX if drifting occurs, perhaps try fiddle with PARTICLE_R where interaction
// is calculated.. it might be that the cells slightly truncate the potential
// number of interactive particles
#define PARTICLE_R (20.0f)
#define PARTICLE_R_SQR (PARTICLE_R*PARTICLE_R)
#define MASS (1.2f)

static float frand(float min, float max)
{
	int r = rand();
	return min + (max - min) * ((float)r / (float)RAND_MAX);
}


static int pot(int value) {
	int pot = 0;
	do { pot++; } while(value >>= 1);
	return pot;
}

static uint32_t murmur32(uint32_t k)
{
	const uint32_t m = 0x5bd1e995;
	const int32_t r = 24;
	uint32_t h = /* seed ^ */ 4;

	k *= m;
	k ^= k >> r;
	k *= m;

	h *= m;
	h ^= k;

	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;

	return h;
}


static uint32_t cell_key(float x, float y)
{
	float s = 1.0f / PARTICLE_R;
	float a = (float)(1<<15);
	float fcx = x * s + a;
	float fcy = y * s + a;
	uint32_t cx = (uint32_t)fcx;
	uint32_t cy = (uint32_t)fcy;
	return (cx&((1<<16)-1)) + (cy << 16);
}

static uint32_t cell_key_hash(uint32_t cell_key)
{
	return murmur32(cell_key);
}


static uint32_t translate_cell_key(uint32_t cell_key, int dx, int dy)
{
	//return cell_key + dx + (dy << 16);
	return (uint32_t)(((int)cell_key) + dx + (dy * (1<<16)));
}

static void particle_update_cell_key(struct particle* p)
{
	p->cell_key = cell_key(p->px, p->py);
}

static void psys_reset_hash(struct psys* ps)
{
	ps->occupied_buckets_count = 0;
	ps->collisions = 0;
	struct particle_hash* ph = &ps->hashes[ps->current_hash_index];
	size_t buckets_bytes = (1<<ps->bucket_count_shift) * sizeof(struct particle_hash_bucket);
	if(ph->buckets == NULL) {
		ph->buckets = malloc(buckets_bytes);
	}
	bzero(ph->buckets, buckets_bytes);
}

static uint32_t mask_bucket_index(struct psys* ps, uint32_t bucket_index)
{
	return bucket_index & ((1<<ps->bucket_count_shift)-1);
}


static void psys_insert_particle(struct psys* ps, struct particle* p)
{
	uint32_t bucket_index = cell_key_hash(p->cell_key);

	bucket_index &= mask_bucket_index(ps, bucket_index);

	struct particle_hash* ph = &ps->hashes[ps->current_hash_index];

	struct particle_hash_bucket* root_bucket = &ph->buckets[bucket_index];
	int insert_offset = root_bucket->insert_offset;

	bucket_index = mask_bucket_index(ps, bucket_index + insert_offset);

	// find empty slot
	while(ph->buckets[bucket_index].particle.cell_key != 0) {
		bucket_index = mask_bucket_index(ps, bucket_index + 1);
		insert_offset++;
		ps->collisions++;
	}

	insert_offset++;

	struct particle* destp = &ph->buckets[bucket_index].particle;
	memcpy(destp, p, sizeof(struct particle));
	particle_update_cell_key(destp);
	destp->density = MASS;

	ps->occupied_buckets[ps->occupied_buckets_count++] = bucket_index;
	root_bucket->insert_offset = insert_offset;
}


static int occupied_buckets_index_compare(const void* va, const void* vb)
{
	int a = *((int*)va);
	int b = *((int*)vb);
	return a - b;
}

static void psys_swap_hashes(struct psys* ps)
{
	ps->current_hash_index ^= 1;
}

void psys_init(struct psys* ps)
{
	bzero(ps, sizeof(struct psys));

	int N = 15000;

	ps->bucket_count_shift = pot(N) + 5; // x32 space per particle - ish
	ps->occupied_buckets_max = N * 2;
	ps->occupied_buckets = malloc(ps->occupied_buckets_max * sizeof(int));

	ps->ppairs = malloc(50 * N * sizeof(struct ppair)); // XXX wishful thinking

	psys_reset_hash(ps);

	struct particle p;
	bzero(&p, sizeof(struct particle));

	for(int i = 0; i < N; i++) {
		p.px = frand(-500, 500);
		p.py = frand(-500, 500);
		//p.vx = frand(-1, 1);
		//p.vy = frand(-1, 1);
		psys_insert_particle(ps, &p);
	}

	// XXX?
	//qsort(ps->occupied_buckets, ps->occupied_buckets_count, sizeof(int), occupied_buckets_index_compare);
}


void psys_step(struct psys* ps)
{
	struct particle_hash* ph = &ps->hashes[ps->current_hash_index];
	int n = ps->occupied_buckets_count;

	ps->ppair_count = 0;

	// neighbour detection and density
	uint32_t d_idx = 0;
	for(int dx = -1; dx <= 1; dx++) {
		for(int dy = -1; dy <= 1; dy++) {
			for(int i = 0; i < n; i++) {
				uint32_t s_idx = ps->occupied_buckets[i];
				struct particle* p = &ph->buckets[s_idx].particle;
				uint32_t neighbour_cell_key = translate_cell_key(p->cell_key, dx, dy);
				d_idx = mask_bucket_index(ps, cell_key_hash(neighbour_cell_key));
				for(;;) {
					struct particle* op = &ph->buckets[d_idx].particle;
					if(op->cell_key == 0) break;
					if(op->cell_key == neighbour_cell_key && s_idx < d_idx) {
						float ddx = op->px - p->px;
						float ddy = op->py - p->py;
						float dsqr = ddx*ddx + ddy*ddy;
						if(dsqr < PARTICLE_R_SQR) {
							// calculate density
							float c = (PARTICLE_R_SQR - dsqr) * (1.0f / PARTICLE_R_SQR);
							float d = c*c*c;
							p->density += d * MASS;
							op->density += d * MASS;

							// add particle pair
							struct ppair* ppair = &ps->ppairs[ps->ppair_count++];
							ppair->i = s_idx;
							ppair->j = d_idx;
						}
					}
					d_idx = mask_bucket_index(ps, d_idx + 1);
				}
			}
		}
	}

	printf("ppair_count: %d\n", ps->ppair_count);

	// calculate forces
	for(int i = 0; i < ps->ppair_count; i++) {
		struct ppair* ppair = &ps->ppairs[i];
		struct particle* particle_i = &ph->buckets[ppair->i].particle;
		struct particle* particle_j = &ph->buckets[ppair->j].particle;
		float dx = particle_i->px - particle_j->px;
		float dy = particle_i->py - particle_j->py;
		float dsqr = dx*dx + dy*dy;
		float r = sqrtf(dsqr);
		float c = 1.0f - r * (1.0f / PARTICLE_R);
		const float d0 = 3.0f; // REST DENSITY
		float pressure_i = particle_i->density - d0;
		float pressure_j = particle_j->density - d0;
		float pterm = (0.5f * c * (pressure_i + pressure_j)) / r;
		float dterm = c / (particle_i->density * particle_j->density);
		//const float step = 0.1f;
		//float pd = pterm * dterm * step;
		//float fx = dx * pd;
		//float fy = dy * pd;

		float visc = 0.05f;
		float fx = (pterm * dx + visc * (particle_j->vxe - particle_i->vxe)) * dterm;
		float fy = (pterm * dy + visc * (particle_j->vye - particle_i->vye)) * dterm;


		// NOTE when pressure is positive, particles are pushed away
		// from each other. opposite holds too; negative pressure
		// attracts particles

		float m1 = 1.0f / MASS;
		fx *= m1;
		fy *= m1;

		// apply forces
		particle_i->vx += fx;
		particle_i->vy += fy;
		particle_j->vx -= fx;
		particle_j->vy -= fy;

		particle_i->vxe += fx;
		particle_i->vye += fy;
		particle_j->vxe -= fx;
		particle_j->vye -= fy;
	}

	// velocity step + rehash
	psys_swap_hashes(ps);
	psys_reset_hash(ps);
	for(int i = 0; i < n; i++) {
		struct particle* p = &ph->buckets[ps->occupied_buckets[i]].particle;

		// gravity
		p->vy += 0.001f;

		// bounds
		if(p->px < -1920/2 || p->px > 1920/2) p->vx = -p->vx * 0.5f;
		if(p->py < -1080/2 || p->py > 1080/2) p->vy = -p->vy * 0.5f;

		p->vxe = (2.0f * p->vx + p->vxe) * 0.5f; // XXX not completely right is it?
		p->vye = (2.0f * p->vy + p->vye) * 0.5f;

		p->px += p->vx;
		p->py += p->vy;
		psys_insert_particle(ps, p);
	}

	// XXX do this? seems like a drop in the water concerning performance.
	// do a cell_key compare too?
	qsort(ps->occupied_buckets, ps->occupied_buckets_count, sizeof(int), occupied_buckets_index_compare);

	printf("collisions: %d\n", ps->collisions);
}

void psys_draw(struct psys* ps)
{
	struct particle_hash* ph = &ps->hashes[ps->current_hash_index];
	glBegin(GL_POINTS);
	for(int i = 0; i < ps->occupied_buckets_count; i++) {
		struct particle* p = &ph->buckets[ps->occupied_buckets[i]].particle;
		glVertex2f(p->px, p->py);
	}
	glEnd();
}


