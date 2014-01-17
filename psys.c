#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>

#include <glew.h>

#include "psys.h"
#include "mud.h"

#define DEG2RAD(x) ((x) * (1.0 / 180.0 * M_PI))

#define SOLID_DIRTY_RGBA (1<<0)
#define SOLID_DIRTY_MASS (1<<1)

#define SOLID_DIRTY_ALL (SOLID_DIRTY_RGBA | SOLID_DIRTY_MASS)


// XXX if drifting occurs, perhaps try fiddle with PARTICLE_R where interaction
// is calculated.. it might be that the cells slightly truncate the potential
// number of interactive particles
#define PARTICLE_R (20.0f)
#define PARTICLE_R_SQR (PARTICLE_R*PARTICLE_R)
#define PARTICLE_MASS (1.2f)



static float dot2d(float ax, float ay, float bx, float by)
{
	return ax*bx + ay*by;
}

/*
static float cross2d(float ax, float ay, float bx, float by)
{
	return ax*by - ay*bx;
}
*/


static uint8_t solid_type_at_point(struct solid* solid, int x, int y)
{
	if(x < 0 || y < 0 || x >= solid->b_width || y >= solid->b_height) {
		return 0;
	} else {
		return *(solid->b_type + x + y * solid->b_width);
	}
}


static float solid_mass_at_point(struct solid* solid, int x, int y)
{
	uint8_t type = solid_type_at_point(solid, x, y);
	switch(type) {
		case 1: return 0.03f;
		default: return 0;
	}
}

static void solid_update_transform(struct solid* solid)
{
	float s = sinf(DEG2RAD(solid->r));
	float c = cosf(DEG2RAD(solid->r));

	float zx = solid->cx * c - solid->cy * s;
	float zy = solid->cx * s + solid->cy * c;

	solid->tx_x0 = solid->px - zx + solid->cx;
	solid->tx_y0 = solid->py - zy + solid->cy;

	//printf("uptx: %f, %f\n", solid->tx_x0, solid->tx_y0);

	solid->tx_u = c;
	solid->tx_v = s;

	// TODO update AABB?
}


/*
static void solid_tx_world_to_local_f(struct solid* solid, float wx, float wy, float* lx, float* ly)
{
	float x0 = wx - solid->tx_x0;
	float y0 = wy - solid->tx_y0;
	*lx = dot2d(solid->tx_u, solid->tx_v, x0, y0);
	*ly = cross2d(solid->tx_u, solid->tx_v, x0, y0);
}

static void solid_tx_world_to_local_i(struct solid* solid, float wx, float wy, int* lx, int* ly)
{
	float lxf;
	float lyf;
	solid_tx_world_to_local_f(solid, wx, wy, &lxf, &lyf);
	// XXX +0.5f?
	// XXX truncate problem around 0
	*lx = (int)lxf;
	*ly = (int)lyf;
}
*/

/*
static int solid_world_point_inside(struct solid* solid, float wx, float wy)
{
	float lx, ly;
	solid_tx_world_to_local_f(solid, wx, wy, &lx, &ly);
	if(lx < 0) return 0;
	if(ly < 0) return 0;
	if(lx >= solid->b_width) return 0;
	if(ly >= solid->b_height) return 0;
	return 1;
}
*/

static int clamp_int(int i, int min, int max)
{
	if(i < min) return min;
	if(i > max) return max;
	return i;
}

static int solid_normal_at_local_point(struct solid* solid, int px, int py, float* nx, float* ny)
{
	const int r = 4;
	int x0 = clamp_int(px - r, 0, solid->b_width - 1);
	int y0 = clamp_int(py - r, 0, solid->b_height - 1);
	int x1 = clamp_int(px + r, 0, solid->b_width - 1);
	int y1 = clamp_int(py + r, 0, solid->b_height - 1);

	int sx = 0;
	int sy = 0;

	uint8_t* t = solid->b_type + x0 + y0 * solid->b_width;
	int j = solid->b_width - (x1 - x0 + 1);
	for(int y = y0; y <= y1; y++) {
		for(int x = x0; x <= x1; x++) {
			uint8_t* t = solid->b_type + x + y * solid->b_width;
			if((*t) != 0) {
				sx += (x-px);
				sy += (y-py);
			}
			t++;
		}
		t += j;
	}

	if(sx == 0 && sy == 0) {
		*nx = 0.0;
		*ny = 0.0;
		return 0;
	} else {
		*nx = -sx;
		*ny = -sy;
		return 1;
	}
}

/*
int solid_normal_at_world_point(struct solid* solid, float px, float py, float* nx, float* ny)
{
	// transform world to local
	float x0 = px - solid->tx_x0;
	float y0 = py - solid->tx_y0;
	float lx = dot2d(solid->tx_u, solid->tx_v, x0, y0);
	float ly = cross2d(solid->tx_u, solid->tx_v, x0, y0);
	int lxi = (int)lx;
	int lyi = (int)ly;

	//solid_draw(solid, lxi, lyi);
	float nnx, nny;

	// find normal
	int res = solid_normal_at_local_point(solid, lxi, lyi, &nnx, &nny);

	// rotate back to world orientation
	(*nx) = nnx * solid->tx_u - nny * solid->tx_v;
	(*ny) = nnx * solid->tx_v + nny * solid->tx_u;

	return res;
}
*/

static void tx_vector(float c, float s, float* u, float* v)
{
	float u0 = *u;
	float v0 = *v;
	float u1 = u0 * c - v0 * s;
	float v1 = u0 * s + v0 * c;
	*u = u1;
	*v = v1;
}

static void solid_tx_vector_local_to_world(struct solid* solid, float* u, float* v)
{
	tx_vector(solid->tx_u, solid->tx_v, u, v);
}

static void solid_tx_vector_world_to_local(struct solid* solid, float* u, float* v)
{
	tx_vector(solid->tx_u, -solid->tx_v, u, v);
}

static int solid_find_escape(struct solid* solid, float x0, float y0, float* x1, float* y1) {
	// search in a spiral until empty spot is found

	int r = 1;

	int cx = x0;
	int cy = y0;

	for(;;) {
		for(int x = (cx-r); x <= (cx+r); x++) {
			if(solid_type_at_point(solid, x, cy-r) == 0) {
				(*x1) = x0 + (x-cx);
				(*y1) = y0 - r;
				return r;
			}
		}
		for(int y = (cy-r); y <= (cy+r); y++) {
			if(solid_type_at_point(solid, cx+r, y) == 0) {
				(*x1) = x0 + r;
				(*y1) = y0 + (y-cy);
				return r;
			}
		}
		for(int x = (cx+r); x >= (cx-r); x--) {
			if(solid_type_at_point(solid, x, cy+r) == 0) {
				(*x1) = x0 + (x-cx);
				(*y1) = y0 + r;
				return r;
			}
		}
		for(int y = (cy+r); y >= (cy-r); y--) {
			if(solid_type_at_point(solid, cx-r, y) == 0) {
				(*x1) = x0 - r;
				(*y1) = y0 + (y-cy);
				return r;
			}
		}
		r++;
	}
}

static void solid_tx_point_world_to_local(struct solid* solid, float* x, float* y)
{
	float lx = (*x) - solid->tx_x0;
	float ly = (*y) - solid->tx_y0;
	*x = dot2d(solid->tx_u, solid->tx_v, lx, ly);
	*y = dot2d(-solid->tx_v, solid->tx_u, lx, ly);
}

static void solid_tx_point_local_to_world(struct solid* solid, float* x, float* y)
{
	solid_tx_vector_local_to_world(solid, x, y);
	(*x) += solid->tx_x0;
	(*y) += solid->tx_y0;
}

static int solid_particle_impulse_response(struct solid* solid, struct particle* particle)
{
	// TODO AABB test

	float lx = particle->px;
	float ly = particle->py;
	solid_tx_point_world_to_local(solid, &lx, &ly);
	int lxi = (int)lx;
	int lyi = (int)ly;

	// check rect
	if(lxi < 0) return 0;
	if(lyi < 0) return 0;
	if(lxi >= solid->b_width) return 0;
	if(lyi >= solid->b_height) return 0;

	// check type cell
	uint8_t t = solid_type_at_point(solid, lxi, lyi);
	if(t == 0) return 0;

	// ok, collision!

	///////////////////////////////////////////////////////
	// calculate impact velocity vector in respect to solid

	// particle velocity in respect to solid
	float pvx = particle->vx;
	float pvy = particle->vy;
	solid_tx_vector_world_to_local(solid, &pvx, &pvy);

	// solid linear velocity in respect to itself (its vx/vy are world
	// velocities)
	float svx = solid->vx;
	float svy = solid->vy;
	solid_tx_vector_world_to_local(solid, &svx, &svy);

	// solid rotational velocity at impact point in respect to itself
	float dx = lx - solid->cx;
	float dy = ly - solid->cy;
	float rvr = DEG2RAD(solid->vr);
	float rvx = -dy * rvr;
	float rvy = dx * rvr;

	// impact velocity in respect to solid
	float ivx = pvx - (svx + rvx);
	float ivy = pvy - (svy + rvy);

	///////////////////////////////////////////////////////

	// estimate impact point _outside_ solid by walking in opposite
	// direction of impact velocity vector
	float px, py;
	solid_find_escape(solid, lx, ly, &px, &py);

	// sample a surface normal
	float nx;
	float ny;
	int found_normal = solid_normal_at_local_point(solid, (int)px, (int)py, &nx, &ny);

	float vn = dot2d(ivx, ivy, nx, ny);

	//printf("DEBUG found_normal = %d at (%d,%d) / vn = %f\n", found_normal, (int)px, (int)py, vn);

	// if a normal was found and impact velocity vector point towards the
	// surface, then we have a "valid impact" and can calculate impulse
	if(found_normal && vn < 0) {
		float nn = nx*nx + ny*ny;
		float m1 = solid->inv_m + (1.0f / PARTICLE_MASS);
		float rn = dot2d(-dy, dx, nx, ny);
		float rnI = rn * rn * solid->inv_I;
		float j = (-2.0f*vn) / (nn * m1 + rnI);

		// apply torque
		solid->fvr -= rn * j; // XXX is this right?

		// linear impulse in respect to world...
		solid_tx_vector_local_to_world(solid, &nx, &ny);
		solid->fx -= nx * j;
		solid->fy -= ny * j;

		particle->vx += nx * j * (1.0 / PARTICLE_MASS);
		particle->vy += ny * j * (1.0 / PARTICLE_MASS);

	}

	solid_tx_point_local_to_world(solid, &px, &py);
	particle->px = px;
	particle->py = py;

	return 1;
}


static void solid_update_dirty(struct solid* solid)
{
	// TODO reupload textures to opengl? (add some dirty-flags maybe)
	// TODO update AABB?

	if(solid->dirty_flags == 0) return;

	if(solid->dirty_flags & SOLID_DIRTY_RGBA) {
		if(solid->gl_rgba != 0) {
		}
		glGenTextures(1, &solid->gl_rgba);
		glBindTexture(GL_TEXTURE_2D, solid->gl_rgba);
		glTexImage2D(
			GL_TEXTURE_2D,
			0, // levels
			4, // bytes per texel
			solid->b_width,
			solid->b_height,
			0,
			GL_RGBA,
			GL_UNSIGNED_BYTE,
			solid->b_rgba);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		solid->dirty_flags &= ~SOLID_DIRTY_RGBA;
	}

	if(solid->dirty_flags & SOLID_DIRTY_MASS) {
		float zmx = 0;
		float zmy = 0;
		float mass = 0;
		for(int y = 0; y < solid->b_height; y++) {
			for(int x = 0; x < solid->b_width; x++) {
				float pm = solid_mass_at_point(solid, x, y);
				mass += pm;
				zmx += pm * ((float)x + 0.5f);
				zmy += pm * ((float)y + 0.5f);
			}
		}
		solid->m = mass;
		solid->inv_m = 1.0f / mass;
		solid->zmx = zmx;
		solid->zmy = zmy;
		solid->cx = zmx / mass;
		solid->cy = zmy / mass;

		float I = 0;
		for(int y = 0; y < solid->b_height; y++) {
			for(int x = 0; x < solid->b_width; x++) {
				float pm = solid_mass_at_point(solid, x, y);
				if(pm > 0.0f) {
					float rx = ((float)x + 0.5f) - solid->cx;
					float ry = ((float)y + 0.5f) - solid->cy;
					float r2 = rx*rx + ry*ry;
					I += pm * r2;
				}
			}
		}
		solid->I = I;
		solid->inv_I = 1.0f / I;

		solid->dirty_flags &= ~SOLID_DIRTY_MASS;
	}
}

/*
struct solid* solid_new(int width, int height)
{
	struct solid* solid = malloc(sizeof(struct solid));
	bzero(solid, sizeof(struct solid));
	solid->b_width = width;
	solid->b_height = height;

	int n = width * height;
	solid->b_rgba = malloc(n * sizeof(uint32_t));
	solid->b_type = malloc(n * sizeof(uint8_t));

	return solid;
}
*/

struct solid* solid_load(const char* id)
{
	struct solid* solid = malloc(sizeof(struct solid));
	bzero(solid, sizeof(struct solid));

	mud_load_png_rgba(id, (void**) &solid->b_rgba, &solid->b_width, &solid->b_height);

	int n = solid->b_width * solid->b_height;

	solid->b_type = (uint8_t*)malloc(n);

	for(int i = 0; i < n; i++) {
		uint32_t* tx = solid->b_rgba + i;
		uint8_t alpha = ((*tx)>>24) & 0xff;
		uint8_t* t = solid->b_type + i;
		(*t) = alpha == 0xff ? 1 : 0;

	}


	// XXX DEBUG
	solid->py = -500;
	solid->r = 90;
	solid->vx = -0.2f;
	solid->vy = 0.13f;
	solid->vr = -0.1f;

	solid->dirty_flags = SOLID_DIRTY_ALL;
	solid_update_dirty(solid);
	solid_update_transform(solid);

	return solid;
}

/*
void solid_destroy(struct solid* solid)
{
	free(solid->b_type);
	free(solid->b_rgba);
	free(solid);
}
*/



static void solid_step(struct solid* solid)
{
	solid->vx += solid->fx * solid->inv_m;
	solid->vy += solid->fy * solid->inv_m;
	solid->fx = 0.0f;
	solid->fy = 0.0f;
	solid->px += solid->vx;
	solid->py += solid->vy;

	solid->vr += solid->fvr * solid->inv_I;
	solid->r += solid->vr;
	solid->fvr = 0.0f;

	solid_update_transform(solid);

	// gravity
	solid->vy += 0.001f;
}


static void solid_draw(struct solid* solid)
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(solid->px + solid->cx, solid->py + solid->cy, 0);
	glRotatef(solid->r, 0, 0, 1);
	glTranslatef(-solid->cx, -solid->cy, 0);

	glColor4f(1,1,1,1);
	glBindTexture(GL_TEXTURE_2D, solid->gl_rgba);
	glBegin(GL_QUADS);
	glTexCoord2f(0, 0);
	glVertex2f(0, 0);
	glTexCoord2f(1, 0);
	glVertex2f(solid->b_width, 0);
	glTexCoord2f(1, 1);
	glVertex2f(solid->b_width, solid->b_height);
	glTexCoord2f(0, 1);
	glVertex2f(0, solid->b_height);
	glEnd();

	// DEBUG: draw transform
	/*
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glColor4f(1,0,1.0,1.0);
	glDisable(GL_TEXTURE_2D);

	float x0 = solid->tx_x0;
	float y0 = solid->tx_y0;

	float ax = solid->tx_u * solid->b_width;
	float ay = solid->tx_v * solid->b_width;

	float bx = -solid->tx_v * solid->b_height;
	float by = solid->tx_u * solid->b_height;

	glBegin(GL_LINE_LOOP);
	glVertex2f(x0, y0);
	glVertex2f(x0 + ax, y0 + ay);
	glVertex2f(x0 + ax + bx, y0 + ay + by);
	glVertex2f(x0 + bx, y0 + by);
	glEnd();
	*/
}




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
	destp->density = PARTICLE_MASS;

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
		p.py = frand(0, 500);
		//p.vx = frand(-1, 1);
		//p.vy = frand(-1, 1);
		psys_insert_particle(ps, &p);
	}

	// XXX?
	//qsort(ps->occupied_buckets, ps->occupied_buckets_count, sizeof(int), occupied_buckets_index_compare);

	ps->solids = solid_load("thing.png");
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
							//float d = c*c; // works too
							//float d = c*c*c*c; // and this
							p->density += d * PARTICLE_MASS;
							op->density += d * PARTICLE_MASS;

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

	//printf("ppair_count: %d\n", ps->ppair_count);

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

		float m1 = 1.0f / PARTICLE_MASS;
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

		p->vxe = (2.0f * p->vx + p->vxe) * 0.5f; // XXX not completely right is it? (but it works!)
		p->vye = (2.0f * p->vy + p->vye) * 0.5f;

		p->px += p->vx;
		p->py += p->vy;

		struct solid* solid = ps->solids;
		while(solid) {
			solid_particle_impulse_response(solid, p);
			solid = solid->next;
		}

		psys_insert_particle(ps, p);
	}

	// XXX do this? seems like a drop in the water concerning performance.
	// do a cell_key compare too?
	qsort(ps->occupied_buckets, ps->occupied_buckets_count, sizeof(int), occupied_buckets_index_compare);

	struct solid* solid = ps->solids;
	while(solid) {
		solid_step(solid);
		solid = solid->next;
	}

	//printf("collisions: %d\n", ps->collisions);
}

void psys_draw(struct psys* ps)
{
	glColor4f(1,1,1,1);
	glDisable(GL_TEXTURE_2D);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	struct particle_hash* ph = &ps->hashes[ps->current_hash_index];
	glBegin(GL_POINTS);
	for(int i = 0; i < ps->occupied_buckets_count; i++) {
		struct particle* p = &ph->buckets[ps->occupied_buckets[i]].particle;
		glVertex2f(p->px, p->py);
	}
	glEnd();

	glEnable(GL_TEXTURE_2D);
	struct solid* solid = ps->solids;
	while(solid) {
		solid_draw(solid);
		solid = solid->next;
	}
}

