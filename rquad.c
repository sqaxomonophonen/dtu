#include <string.h>
#include <math.h>

#include "rquad.h"
#include "a.h"

struct rquad_split {
	float x0;
	float y0;
	float u;
	float v;
};

static float split_dist(struct rquad_split* split, struct rquad_point* point)
{
	float nx = -split->v;
	float ny = split->u;
	float dx = point->x - split->x0;
	float dy = point->y - split->y0;
	return dx*nx + dy*ny;
}

static float lerp(float a, float b, float t)
{
	return a + (b-a) * t;
}

static float cross(float ax, float ay, float bx, float by)
{
	return ax * by - ay * bx;
}

static void clip(
	struct rquad_point* p,
	struct rquad_point* p0,
	struct rquad_point* p1,
	struct rquad_split* s
) {
	float denom = cross(p1->x - p0->x, p1->y - p0->y, s->u, s->v);
	float numer = cross(s->x0 - p0->x, s->y0 - p0->y, s->u, s->v);
	float t = numer/denom;
	p->x = lerp(p0->x, p1->x, t);
	p->y = lerp(p0->y, p1->y, t);
	p->u = lerp(p0->u, p1->u, t);
	p->v = lerp(p0->v, p1->v, t);
}

int rquad_init(struct rquad* rquad) {
	bzero(&rquad->_internal, sizeof(struct rquad_internal));

	const float x0 = rquad->x0;
	const float y0 = rquad->y0;
	const float xdw = rquad->u * rquad->src_width;
	const float ydw = rquad->v * rquad->src_width;
	const float xdh = -rquad->v * rquad->src_height;
	const float ydh = rquad->u * rquad->src_height;
	const float ESx0 = 0.0f;
	const float ESy0 = 0.0f;
	const float ESx1 = 1.0f - ESx0;
	const float ESy1 = 1.0f - ESy0;

	// with EPSILON -- XXX is it required?
	/*
	const float E = 0.1f;
	const float E2 = E*2.0f;
	const float x0 = rquad->x0 + E;
	const float y0 = rquad->y0 + E;
	const float xdw = rquad->u * (rquad->src_width - E2);
	const float ydw = rquad->v * (rquad->src_width - E2);
	const float xdh = -rquad->v * (rquad->src_height - E2);
	const float ydh = rquad->u * (rquad->src_height - E2);
	const float ESx0 = E / rquad->src_width;
	const float ESy0 = E / rquad->src_height;
	const float ESx1 = 1.0f - ESx0;
	const float ESy1 = 1.0f - ESy0;
	*/

	struct rquad_point points[RQUAD_MAX_POINTS] = {
		{x0, y0, ESx0, ESy0},
		{x0 + xdw, y0 + ydw, ESx1, ESy0},
		{x0 + xdw + xdh, y0 + ydw + ydh, ESx1, ESy1},
		{x0 + xdh, y0 + ydh, ESx0, ESy1},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0}
	};

	int point_count = 4;

	struct rquad_split splits[4] = {
		{0, 0, 1, 0},
		{0, 0, 0, -1},
		{rquad->dst_width, 0, 0, 1},
		{0, rquad->dst_height, -1, 0}
	};

	// clip sides to bounds
	for(int split_index = 0; split_index < 4; split_index++) {
		struct rquad_split* split = &splits[split_index];
		int side[RQUAD_MAX_POINTS];
		for(int p = 0; p < point_count; p++) {
			struct rquad_point* point = &points[p];
			float d = split_dist(split, point);
			side[p] = d >= 0.0f;
		}

		struct rquad_point new_points[RQUAD_MAX_POINTS];
		int new_point_count = 0;

		for(int i0 = 0; i0 < point_count; i0++) {
			int i1 = _rquad_wrap_index(i0 + 1, point_count);
			struct rquad_point* p0 = &points[i0];
			struct rquad_point* p1 = &points[i1];
			int side0 = side[i0];
			int side1 = side[i1];

			size_t psz = sizeof(struct rquad_point);
			if(!side0 && !side1) {
				continue;
			} else if(side0 && side1) {
				memcpy(&new_points[new_point_count++], p0, psz);
			} else if(side0 && !side1) {
				memcpy(&new_points[new_point_count++], p0, psz);
				struct rquad_point p;
				clip(&p, p0, p1, split);
				memcpy(&new_points[new_point_count++], &p, psz);
			} else if(!side0 && side1) {
				struct rquad_point p;
				clip(&p, p0, p1, split);
				memcpy(&new_points[new_point_count++], &p, psz);
			}
		}

		if(new_point_count == 0) {
			// all points culled; no point in continuing
			return 0;
		} else {
			memcpy(
				points,
				new_points,
				sizeof(struct rquad_point) * new_point_count);
			point_count = new_point_count;
		}
	}

	if(point_count < 3) {
		// degenerate case
		return 0;
	}

	// prep rquad
	struct rquad_internal* rqi = &rquad->_internal;
	rqi->point_count = point_count;
	memcpy(
		rqi->points,
		points,
		sizeof(struct rquad_point) * point_count);
	rqi->triangles_left = point_count - 2;

	return 1;
}

void _rquad_setup_triangle(struct rquad* rq)
{
	struct rquad_internal* rqi = &rq->_internal;

	rqi->sub_triangles_left = 0;

	int triangle_index = rqi->point_count - 2 - rqi->triangles_left;

	struct rquad_point* v0 = &rqi->points[0];
	struct rquad_point* v1 = &rqi->points[_rquad_wrap_index(
		1+triangle_index,
		rqi->point_count)];
	struct rquad_point* v2 = &rqi->points[_rquad_wrap_index(
		2+triangle_index,
		rqi->point_count)];

	const int snap_mask = ~((_RQUAD_FIXED_ONE / (1 << _RQUAD_SUB_PIXEL_BITS)) - 1);

	int fy0 = _rquad_float_to_fixed(v0->y - 0.5f) & snap_mask;
	int fy1 = _rquad_float_to_fixed(v1->y - 0.5f) & snap_mask;
	int fy2 = _rquad_float_to_fixed(v2->y - 0.5f) & snap_mask;

	struct rquad_point* v_min;
	struct rquad_point* v_mid;
	struct rquad_point* v_max;
	int v_min_fy, v_mid_fy, v_max_fy;

	if(fy0 <= fy1) {
		if(fy1 <= fy2) {
			/* y0 <= y1 <= y2 */
			v_min = v0; v_mid = v1; v_max = v2;
			v_min_fy = fy0; v_mid_fy = fy1; v_max_fy = fy2;
		} else if(fy2 <= fy0) {
			/* y2 <= y0 <= y1 */
			v_min = v2; v_mid = v0; v_max = v1;
			v_min_fy = fy2; v_mid_fy = fy0; v_max_fy = fy1;
		} else {
			/* y0 <= y2 <= y1 */
			v_min = v0; v_mid = v2; v_max = v1;
			v_min_fy = fy0; v_mid_fy = fy2; v_max_fy = fy1;
		}
	} else {
		if(fy0 <= fy2) {
			/* y1 <= y0 <= y2 */
			v_min = v1; v_mid = v0; v_max = v2;
			v_min_fy = fy1; v_mid_fy = fy0; v_max_fy = fy2;
		} else if(fy2 <= fy1) {
			/* y2 <= y1 <= y0 */
			v_min = v2; v_mid = v1; v_max = v0;
			v_min_fy = fy2; v_mid_fy = fy1; v_max_fy = fy0;
		} else {
			/* y1 <= y2 <= y0 */
			v_min = v1; v_mid = v2; v_max = v0;
			v_min_fy = fy1; v_mid_fy = fy2; v_max_fy = fy0;
		}
	}

	int v_min_fx = _rquad_float_to_fixed(v_min->x + 0.5f) & snap_mask;
	int v_mid_fx = _rquad_float_to_fixed(v_mid->x + 0.5f) & snap_mask;
	int v_max_fx = _rquad_float_to_fixed(v_max->x + 0.5f) & snap_mask;

	rqi->e_maj.v0 = v_min; rqi->e_maj.v1 = v_max;
	rqi->e_top.v0 = v_mid; rqi->e_top.v1 = v_max;
	rqi->e_bot.v0 = v_min; rqi->e_bot.v1 = v_mid;

	rqi->e_maj.dx = _rquad_fixed_to_float(v_max_fx - v_min_fx);
	rqi->e_maj.dy = _rquad_fixed_to_float(v_max_fy - v_min_fy);

	rqi->e_top.dx = _rquad_fixed_to_float(v_max_fx - v_mid_fx);
	rqi->e_top.dy = _rquad_fixed_to_float(v_max_fy - v_mid_fy);

	rqi->e_bot.dx = _rquad_fixed_to_float(v_mid_fx - v_min_fx);
	rqi->e_bot.dy = _rquad_fixed_to_float(v_mid_fy - v_min_fy);

	float area = rqi->e_maj.dx * rqi->e_bot.dy - rqi->e_bot.dx * rqi->e_maj.dy;

	if(!isfinite(area) || area == 0.0f) {
		// ignore triangle, and skip to next
		return;
	}

	float one_over_area = 1.0f / area;

	rqi->e_maj.fsy = _rquad_fixed_ceil(v_min_fy);
	rqi->e_maj.row_count = _rquad_fixed_to_int(_rquad_fixed_ceil(v_max_fy - rqi->e_maj.fsy));
	if(rqi->e_maj.row_count > 0) {
		rqi->e_maj.dxdy = rqi->e_maj.dx / rqi->e_maj.dy;
		rqi->e_maj.fdxdy = _rquad_float_to_fixed(rqi->e_maj.dxdy);
		rqi->e_maj.adjy = (float) (rqi->e_maj.fsy - v_min_fy);  /* SCALED! */
		rqi->e_maj.fx0 = v_min_fx;
		rqi->e_maj.fsx = rqi->e_maj.fx0 + (int) (rqi->e_maj.adjy * rqi->e_maj.dxdy);
	} else {
		return;
	}

	rqi->e_top.fsy = _rquad_fixed_ceil(v_mid_fy);
	rqi->e_top.row_count = _rquad_fixed_to_int(_rquad_fixed_ceil(v_max_fy - rqi->e_top.fsy));
	if (rqi->e_top.row_count > 0) {
		rqi->e_top.dxdy = rqi->e_top.dx / rqi->e_top.dy;
		rqi->e_top.fdxdy = _rquad_float_to_fixed(rqi->e_top.dxdy);
		rqi->e_top.adjy = (float) (rqi->e_top.fsy - v_mid_fy); /* SCALED! */
		rqi->e_top.fx0 = v_mid_fx;
		rqi->e_top.fsx = rqi->e_top.fx0 + (int) (rqi->e_top.adjy * rqi->e_top.dxdy);
	}

	rqi->e_bot.fsy = _rquad_fixed_ceil(v_min_fy);
	rqi->e_bot.row_count = _rquad_fixed_to_int(_rquad_fixed_ceil(v_mid_fy - rqi->e_bot.fsy));
	if (rqi->e_bot.row_count > 0) {
		rqi->e_bot.dxdy = rqi->e_bot.dx / rqi->e_bot.dy;
		rqi->e_bot.fdxdy = _rquad_float_to_fixed(rqi->e_bot.dxdy);
		rqi->e_bot.adjy = (float) (rqi->e_bot.fsy - v_min_fy);  /* SCALED! */
		rqi->e_bot.fx0 = v_min_fx;
		rqi->e_bot.fsx = rqi->e_bot.fx0 + (int) (rqi->e_bot.adjy * rqi->e_bot.dxdy);
	}

	rqi->scan_from_left_to_right = one_over_area < 0.0f;

	float e_maj_ds = (v_max->u - v_min->u) * rq->src_width;
	float e_bot_ds = (v_mid->u - v_min->u) * rq->src_width;
	float e_maj_dt = (v_max->v - v_min->v) * rq->src_height;
	float e_bot_dt = (v_mid->v - v_min->v) * rq->src_height;

	rqi->sx_step = one_over_area * (e_maj_ds * rqi->e_bot.dy - rqi->e_maj.dy * e_bot_ds);
	rqi->sy_step = one_over_area * (rqi->e_maj.dx * e_bot_ds - e_maj_ds * rqi->e_bot.dx);
	rqi->tx_step = one_over_area * (e_maj_dt * rqi->e_bot.dy - rqi->e_maj.dy * e_bot_dt);
	rqi->ty_step = one_over_area * (rqi->e_maj.dx * e_bot_dt - e_maj_dt * rqi->e_bot.dx);
	rqi->fsx_step = _rquad_float_to_fixed(rqi->sx_step);
	rqi->ftx_step = _rquad_float_to_fixed(rqi->tx_step);

	rqi->sub_triangles_left = 2;
}

void _rquad_setup_sub_triangle(struct rquad* rq)
{
	struct rquad_internal* rqi = &rq->_internal;

	rqi->rows_left = 0;

	int sub_triangle_index = 2 - rqi->sub_triangles_left;
	ASSERT(sub_triangle_index >= 0 && sub_triangle_index <= 1);

	struct rquad_edge* e_left;
	struct rquad_edge* e_right;
	int setup_left;
	int setup_right;

	if(sub_triangle_index == 0) {
		/* bottom half */
		if(rqi->scan_from_left_to_right) {
			e_left = &rqi->e_maj;
			e_right = &rqi->e_bot;
			rqi->rows_left = e_right->row_count;
			setup_left = 1;
			setup_right = 1;
		} else {
			e_left = &rqi->e_bot;
			e_right = &rqi->e_maj;
			rqi->rows_left = e_left->row_count;
			setup_left = 1;
			setup_right = 1;
		}
	} else {
		/* top half */
		if(rqi->scan_from_left_to_right) {
			e_left = &rqi->e_maj;
			e_right = &rqi->e_top;
			rqi->rows_left = e_right->row_count;
			setup_left = 0;
			setup_right = 1;
		} else {
			e_left = &rqi->e_top;
			e_right = &rqi->e_maj;
			rqi->rows_left = e_left->row_count;
			setup_left = 1;
			setup_right = 0;
		}
		if(!rqi->rows_left) return;
	}

	if(setup_left && e_left->row_count > 0) {
		const struct rquad_point* v_lower = e_left->v0;
		const int fsy = e_left->fsy;
		const int fsx = e_left->fsx;  /* no fractional part */
		const int fx = _rquad_fixed_ceil(fsx);  /* no fractional part */
		const int adjx = (int) (fx - e_left->fx0); /* SCALED! */
		const int adjy = (int) e_left->adjy;      /* SCALED! */

		rqi->f_error = fx - fsx - _RQUAD_FIXED_ONE;
		rqi->fx_left_edge = fsx - _RQUAD_FIXED_EPSILON;
		rqi->fdx_left_edge = e_left->fdxdy;
		int fdx_outer = _rquad_fixed_floor(rqi->fdx_left_edge - _RQUAD_FIXED_EPSILON);
		rqi->fd_error = fdx_outer - rqi->fdx_left_edge + _RQUAD_FIXED_ONE;
		int idx_outer = _rquad_fixed_to_int(fdx_outer);
		float dx_outer = idx_outer;
		rqi->y = _rquad_fixed_to_int(fsy);

		float s0 = v_lower->u * rq->src_width;
		rqi->s_left = (int)(
			s0 * _RQUAD_FIXED_SCALE
			+ rqi->sx_step * adjx
			+ rqi->sy_step * adjy
			) + _RQUAD_FIXED_HALF;
		rqi->ds_outer = _rquad_float_to_fixed(
			rqi->sy_step
			+ dx_outer * rqi->sx_step);
		float t0 = v_lower->v * rq->src_height;
		rqi->t_left = (int)(
			t0 * _RQUAD_FIXED_SCALE
			+ rqi->tx_step * adjx
			+ rqi->ty_step * adjy
			) + _RQUAD_FIXED_HALF;
		rqi->dt_outer = _rquad_float_to_fixed(
			rqi->ty_step
			+ dx_outer * rqi->tx_step);

	}

	if(setup_right && e_right->row_count > 0) {
		rqi->fx_right_edge = e_right->fsx - _RQUAD_FIXED_EPSILON;
		rqi->fdx_right_edge = e_right->fdxdy;
	}

	if(!rqi->rows_left) return;

	rqi->ds_inner = rqi->ds_outer + rqi->fsx_step;
	rqi->dt_inner = rqi->dt_outer + rqi->ftx_step;
}


#ifdef RQUAD_UNIT_TEST
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define T(FN) do { if(FN()) { fprintf(stdout, "%s OK\n", #FN); } else { fprintf(stderr, "%s FAILED!\n", #FN); exit(EXIT_FAILURE); } } while(0);

static int test_blit()
{
	const int w = 11;
	const int h = 7;
	const int N = w * h;

	struct rquad rq;
	rq.dst_width = w;
	rq.dst_height = h;
	rq.src_width = w;
	rq.src_height = h;
	rq.x0 = 0;
	rq.y0 = 0;
	rq.u = 1;
	rq.v = 0;

	int count[N];
	bzero(count, sizeof(count));

	int noni = 0;

	if(rquad_init(&rq)) {
		if(rq._internal.point_count != 4) {
			fprintf(stderr, "point_count != 4 (was %d)\n", rq._internal.point_count);
			return 0;
		}
		while(rquad_step(&rq)) {
			int x = rquad_x(&rq);
			int y = rquad_y(&rq);
			int s = rquad_s(&rq);
			int t = rquad_t(&rq);

			if(x < 0 || x >= w) {
				fprintf(stderr, "x out of bounds! (x=%d)\n", x);
				return 0;
			}
			if(y < 0 || y >= h) {
				fprintf(stderr, "y out of bounds! (y=%d)\n", y);
				return 0;
			}

			if(s < 0 || s >= w) {
				fprintf(stderr, "s out of bounds! (s=%d)\n", s);
				return 0;
			}

			if(t < 0 || t >= h) {
				fprintf(stderr, "t out of bounds! (t=%d)\n", t);
				return 0;
			}

			if(x != s || y != t) {
				fprintf(stderr, "(x,y) vs (s,t) mismatch: (%d, %d) != (%d, %d)\n", x, y, s, t);
				noni++;
			}

			count[x + y * w]++;
		}

		if(noni) {
			return 0;
		}

		int z = 0;
		int o = 0;
		for(int i = 0; i < N; i++) {
			if(count[i] == 0) {
				z++;
			} else if(count[i] != 1) {
				o++;
			}
		}
		if(z != 0 || o != 0) {
			fprintf(stderr, "expected all pixels to be written exactly once, zero times: %d, more than once: %d\n", z, o);
			return 0;
		}
	} else {
		fprintf(stderr, "expected _something_ to rasterize\n");
		return 0;
	}
	return 1;
}

static int test_rotation_blit_90()
{
	const int S = 99;
	const int N = S*S;

	int count[N];

	for(int ri = 0; ri < 4; ri++) {
		int noni = 0;
		bzero(count, sizeof(count));

		struct rquad rq;
		rq.dst_width = S;
		rq.dst_height = S;
		rq.src_width = S;
		rq.src_height = S;
		switch(ri) {
		case 0:
			rq.x0 = 0;
			rq.y0 = 0;
			rq.u = 1;
			rq.v = 0;
			break;
		case 1:
			rq.x0 = S;
			rq.y0 = 0;
			rq.u = 0;
			rq.v = 1;
			break;
		case 2:
			rq.x0 = S;
			rq.y0 = S;
			rq.u = -1;
			rq.v = 0;
			break;
		case 3:
			rq.x0 = 0;
			rq.y0 = S;
			rq.u = 0;
			rq.v = -1;
			break;
		}

		if(!rquad_init(&rq)) {
			fprintf(stderr, "expected to rasterize _something_ (ri=%d)\n", ri);
			return 0;
		}
		while(rquad_step(&rq)) {
			int x = rquad_x(&rq);
			int y = rquad_y(&rq);
			int s = rquad_s(&rq);
			int t = rquad_t(&rq);

			if(x < 0 || x >= S) {
				fprintf(stderr, "(ri=%d) x out of bounds! (x=%d)\n", ri, x);
				return 0;
			}
			if(y < 0 || y >= S) {
				fprintf(stderr, "(ri=%d) y out of bounds! (y=%d)\n", ri, y);
				return 0;
			}

			if(s < 0 || s >= S) {
				fprintf(stderr, "(ri=%d) s out of bounds! (s=%d)\n", ri, s);
				return 0;
			}

			if(t < 0 || t >= S) {
				fprintf(stderr, "(ri=%d) t out of bounds! (t=%d)\n", ri, t);
				return 0;
			}

			int xs;
			int xt;
			switch(ri) {
			case 0:
				xs = x;
				xt = y;
				break;
			case 1:
				xs = y;
				xt = S - x - 1;
				break;
			case 2:
				xs = S - x - 1;
				xt = S - y - 1;
				break;
			case 3:
				xs = S - y - 1;
				xt = x;
				break;

			}

			if(xs != s || xt != t) {
				fprintf(stderr, "(ri=%d) expected (s,t) to be (%d, %d), got (%d, %d)\n", ri, xs, xt, s, t);
				noni++;
			}

			count[x + y * S]++;
		}

		if(noni) {
			return 0;
		}

		int z = 0;
		int o = 0;
		for(int i = 0; i < N; i++) {
			if(count[i] == 0) {
				z++;
			} else if(count[i] != 1) {
				o++;
			}
		}
		if(z != 0 || o != 0) {
			fprintf(stderr, "(ri=%d) expected all pixels to be written exactly once, zero times: %d, more than once: %d\n", ri, z, o);
			return 0;
		}

	}
	return 1;
}

static float randf(float min, float max)
{
	int r = rand();
	float rf = (float)r / (float)RAND_MAX;
	return min + (max-min) * rf;
}

static void dump_rq(struct rquad* rq)
{
	fprintf(stderr, "RQUAD DUMP:\n");
	fprintf(stderr, "  rq.dst_width = %d\n", rq->dst_width);
	fprintf(stderr, "  rq.dst_height = %d\n", rq->dst_height);
	fprintf(stderr, "  rq.src_width = %d\n", rq->src_width);
	fprintf(stderr, "  rq.src_height = %d\n", rq->src_height);
	fprintf(stderr, "  rq.x0 = %f\n", rq->x0);
	fprintf(stderr, "  rq.y0 = %f\n", rq->y0);
	fprintf(stderr, "  rq.u = %f\n", rq->u);
	fprintf(stderr, "  rq.v = %f\n", rq->v);
}

static int test_against_s_t_overflow()
{
	// XXX this is currently kind of a non-test since s/t are being clamped :-/

	const int S = 103;

	int ok = 0;
	int fail = 0;

	for(int i = 0; i < 5000; i++) {
		struct rquad rq;
		rq.dst_width = S;
		rq.dst_height = S;
		rq.src_width = S/2;
		rq.src_height = S/2;
		rq.x0 = ((float)S / 4.0f) + randf(-0.6f, 0.6f);
		rq.y0 = ((float)S / 4.0f) + randf(-0.6f, 0.6f);
		float rot = randf(-0.1f, 0.1f);
		rq.u = cosf(rot);
		rq.v = sinf(rot);

		if(!rquad_init(&rq)) {
			fprintf(stderr, "expected to rasterize _something_\n");
			return 0;
		}

		int bounds_fail = 0;

		while(rquad_step(&rq)) {
			int s = rquad_s(&rq);
			int t = rquad_t(&rq);
			if(s < 0 || s >= rq.src_width || t < 0 || t >= rq.src_height) {
				fprintf(stderr, "(s,t) out of bounds: (%d,%d)\n", s, t);
				bounds_fail++;
			}
		}

		if(bounds_fail) {
			dump_rq(&rq);
			fail++;
		} else {
			ok++;
		}
	}

	if(fail) {
		fprintf(stderr, "%d boundary overflows, %d ok\n", fail, ok);
		return 0;
	}
	return 1;
}

static int test_against_x_y_overflow()
{
	const int S = 101;

	int ok = 0;
	int fail = 0;

	for(int i = 0; i < 5000; i++) {
		struct rquad rq;
		rq.dst_width = S;
		rq.dst_height = S;
		rq.src_width = S*0.9;
		rq.src_height = S*0.9;
		rq.x0 = randf(-3.0f, 3.0f);
		rq.y0 = randf(-3.0f, 3.0f);
		float rot = randf(-1.0f, 1.0f);
		rq.u = cosf(rot);
		rq.v = sinf(rot);

		if(!rquad_init(&rq)) {
			fprintf(stderr, "expected to rasterize _something_\n");
			return 0;
		}

		int bounds_fail = 0;

		while(rquad_step(&rq)) {
			int x = rquad_x(&rq);
			int y = rquad_y(&rq);
			if(x < 0 || x >= rq.dst_width || y < 0 || y >= rq.dst_height) {
				fprintf(stderr, "(x,y) out of bounds: (%d,%d)\n", x, y);
				bounds_fail++;
			}
		}

		if(bounds_fail) {
			dump_rq(&rq);
			fail++;
		} else {
			ok++;
		}
	}

	if(fail) {
		fprintf(stderr, "%d boundary overflows, %d ok\n", fail, ok);
		return 0;
	}
	return 1;
}

static int test_against_overdraw()
{
	const int rows = 7;
	const int columns = 12;

	const int DS = 15;
	const int DN = DS*DS;

	int ok = 0;
	int overdraw = 0;

	for(int d = -400; d <= 400; d++) {
		struct rquad rq;
		rq.dst_width = DS;
		rq.dst_height = DS;
		rq.src_width = columns;
		rq.src_height = rows;
		float r = (float)d / 180.0f * M_PI;
		float c = cosf(r);
		float s = sinf(r);
		float x0 = (DS-rq.src_width)/2;
		float y0 = (DS-rq.src_height)/2;
		rq.x0 = c * (x0 - DS/2) - s * (y0 - DS/2) + DS/2;
		rq.y0 = s * (x0 - DS/2) + c * (y0 - DS/2) + DS/2;
		rq.u = c;
		rq.v = s;

		int buf[DN];
		bzero(buf, sizeof(buf));

		rquad_init(&rq);
		while(rquad_step(&rq)) {
			int x = rquad_x(&rq);
			int y = rquad_y(&rq);
			int s = rquad_s(&rq);
			int t = rquad_t(&rq);

			if(x < 0 || x >= rq.dst_width) {
				fprintf(stderr, "x out of bounds! (x=%d)\n", x);
				return 0;
			}
			if(y < 0 || y >= rq.dst_height) {
				fprintf(stderr, "y out of bounds! (y=%d)\n", y);
				return 0;
			}

			if(s < 0 || s >= rq.src_width) {
				fprintf(stderr, "s out of bounds! (s=%d)\n", s);
				return 0;
			}

			if(t < 0 || t >= rq.src_height) {
				fprintf(stderr, "t out of bounds! (t=%d)\n", t);
				return 0;
			}

			int dsti = x + y * DS;
			buf[dsti]++;
		}

		for(int i = 0; i < DN; i++) {
			if(buf[i] > 1) {
				overdraw++;
			} else {
				ok++;
			}
		}
	}
	if(overdraw) {
		fprintf(stderr, "%d pixels overdrawn, %d ok\n", overdraw, ok);
		return 0;
	}
	return 1;
}

int main(int argc, char** argv)
{
	printf("sizeof(struct rquad) = %zu\n", sizeof(struct rquad));

	srand(time(NULL));

	T(test_blit);
	T(test_rotation_blit_90);
	T(test_against_s_t_overflow);
	T(test_against_x_y_overflow);
	T(test_against_overdraw);

	return EXIT_SUCCESS;
}

#endif // RQUAD_UNIT_TEST


#ifdef RQUAD_VIS_TEST
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/select.h>

static void rotate_a_P()
{
	// XXX really looks like it makes a lot of gaps..

	const char* TEX =
	//	 123456789abc
		"!..######;.." // 1
		"...##;..##;." // 2
		"...##;..##;." // 3
		"...######;.." // 4
		"...##;......" // 5
		"...##;......" // 6
		"...##;.....?";// 7

	const char* TEX2 =
	//	 123456789abc
		"++++++++++++" // 1
		"+**********+" // 2
		"+*========*+" // 3
		"+*########*+" // 4
		"+*========*+" // 5
		"+**********+" // 6
		"++++++++++++";// 7
	const int rows = 7;
	const int columns = 12;

	const int DS = 10;
	const int DN = DS*DS;
	char scr[DN];
	char scr2[DN];

	for(int d = -67; d <= 400; d++) {
		struct rquad rq;
		rq.dst_width = DS;
		rq.dst_height = DS;
		rq.src_width = columns;
		rq.src_height = rows;
		float r = (float)d / 180.0f * M_PI;
		float c = cosf(r);
		float s = sinf(r);
		float x0 = (DS-rq.src_width)/2;
		float y0 = (DS-rq.src_height)/2;
		rq.x0 = c * (x0 - DS/2) - s * (y0 - DS/2) + DS/2;
		rq.y0 = s * (x0 - DS/2) + c * (y0 - DS/2) + DS/2;
		rq.u = c;
		rq.v = s;

		bzero(scr, DN);
		bzero(scr2, DN);

		rquad_init(&rq);
		while(rquad_step(&rq)) {
			int x = rquad_x(&rq);
			int y = rquad_y(&rq);
			int s = rquad_s(&rq);
			int t = rquad_t(&rq);

			if(x < 0 || x >= rq.dst_width) {
				fprintf(stderr, "x out of bounds! (x=%d)\n", x);
				exit(EXIT_FAILURE);
			}
			if(y < 0 || y >= rq.dst_height) {
				fprintf(stderr, "y out of bounds! (y=%d)\n", y);
				exit(EXIT_FAILURE);
			}

			if(s < 0 || s >= rq.src_width) {
				fprintf(stderr, "s out of bounds! (s=%d)\n", s);
				exit(EXIT_FAILURE);
			}

			if(t < 0 || t >= rq.src_height) {
				fprintf(stderr, "t out of bounds! (t=%d)\n", t);
				exit(EXIT_FAILURE);
			}

			int srci = s + t * columns;
			int dsti = x + y * DS;
			scr[dsti] = TEX[srci];
			scr2[dsti] = TEX2[srci];
		}

		//printf("\033[2J\033[H");

		for(int i = 0; i < 2; i++) {
			char* buf;
			if(i == 0) buf = scr; else buf = scr2;
			printf("+");
			for(int i = 0; i < DS; i++) {
				printf("-");
			}
			printf("+\n");

			for(int y = 0; y < DS; y++) {
				printf("|");
				for(int x = 0; x < DS; x++) {
					char p = buf[x + y * DS];
					if(p == 0) {
						printf(" ");
					} else {
						printf("%c", p);
					}
				}
				printf("|\n");
			}

			printf("+");
			for(int i = 0; i < DS; i++) {
				printf("-");
			}
			printf("+\n\n");
		}
		printf("%d degrees\n", d);

		struct timeval tv;
		tv.tv_sec = 0;
		tv.tv_usec = 1000000/40;
		select(0, NULL, NULL, NULL, &tv);
	}
}

int main(int argc, char** argv)
{
	rotate_a_P();
	return EXIT_SUCCESS;
}

#endif // RQUAD_VIS_TEST
