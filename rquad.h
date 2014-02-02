#ifndef __RQUAD_H__
#define __RQUAD_H__

#define RQUAD_MAX_POINTS (8)

struct rquad_point {
	float x;
	float y;
	float u;
	float v;
};

struct rquad_edge {
	struct rquad_point* v0;
	struct rquad_point* v1;
	float dx;
	float dy;
	float dxdy;
	int fdxdy;
	float adjy;
	int fsx;
	int fsy;
	int fx0;
	int row_count;
};

struct rquad_internal {
	int point_count;
	struct rquad_point points[RQUAD_MAX_POINTS];

	// s/t (current texture/src position)
	int fs, ft;
	int fsx_step, ftx_step;
	float sx_step, tx_step, sy_step, ty_step;

	// counters
	int triangles_left;
	int sub_triangles_left;
	int columns_left;
	int rows_left;

	struct rquad_edge e_maj, e_top, e_bot;
	int scan_from_left_to_right;
	int f_error, fd_error;
	int fx_left_edge;
	int fdx_left_edge;
	int fx_right_edge;
	int fdx_right_edge;
	int x, y;
	int s_left, ds_outer, t_left, dt_outer;
	int ds_inner, dt_inner;
};

struct rquad {
	int dst_width;
	int dst_height;

	int src_width;
	int src_height;

	float x0;
	float y0;
	float u;
	float v;

	struct rquad_internal _internal;
};

int rquad_init(struct rquad* rquad);

#define _RQUAD_SUB_PIXEL_BITS 4

static inline int _rquad_iround(float f)
{
	return (int) ((f >= 0.0f) ? (f + 0.5f) : (f - 0.5f));
}

#define _RQUAD_FIXED_FRAC_BITS 11
#define _RQUAD_FIXED_SHIFT     _RQUAD_FIXED_FRAC_BITS
#define _RQUAD_FIXED_ONE       (1 << _RQUAD_FIXED_SHIFT)
#define _RQUAD_FIXED_HALF      (1 << (_RQUAD_FIXED_SHIFT-1))
#define _RQUAD_FIXED_FRAC_MASK (_RQUAD_FIXED_ONE - 1)
#define _RQUAD_FIXED_INT_MASK  (~_RQUAD_FIXED_FRAC_MASK)
#define _RQUAD_FIXED_EPSILON   1
#define _RQUAD_FIXED_SCALE     ((float) _RQUAD_FIXED_ONE)

static inline int _rquad_float_to_fixed(float f)
{
	return _rquad_iround(f * _RQUAD_FIXED_SCALE);
}

static inline int _rquad_fixed_to_int(int i)
{
	return i >> _RQUAD_FIXED_SHIFT;
}

static inline int _rquad_fixed_ceil(int i)
{
	return (i + _RQUAD_FIXED_ONE - _RQUAD_FIXED_EPSILON) & _RQUAD_FIXED_INT_MASK;
}

static inline int _rquad_fixed_floor(int i)
{
	return i & _RQUAD_FIXED_INT_MASK;
}

static inline float _rquad_fixed_to_float(float f)
{
	return f * (1.0f / _RQUAD_FIXED_SCALE);
}

static inline int _rquad_wrap_index(int i, int n)
{
	while(i >= n) i -= n;
	return i;
}

void _rquad_setup_triangle(struct rquad* rq);
void _rquad_setup_sub_triangle(struct rquad* rq);

static inline int rquad_step(struct rquad* rq)
{
	struct rquad_internal* rqi = &rq->_internal;

	//do_col_adv:
	if(rqi->columns_left) {
		rqi->x++;
		rqi->fs += rqi->fsx_step;
		rqi->ft += rqi->ftx_step;
	}

	do_col:
	if(rqi->columns_left) {
		rqi->columns_left--;
		return 1;
	}

	//do_row_adv:
	if(rqi->rows_left) {
		rqi->y++;
		rqi->fx_left_edge += rqi->fdx_left_edge;
		rqi->fx_right_edge += rqi->fdx_right_edge;
		rqi->f_error += rqi->fd_error;
		if(rqi->f_error >= 0) {
			rqi->f_error -= _RQUAD_FIXED_ONE;
			rqi->s_left += rqi->ds_outer;
			rqi->t_left += rqi->dt_outer;
		} else {
			rqi->s_left += rqi->ds_inner;
			rqi->t_left += rqi->dt_inner;
		}
	}

	do_row:
	if(rqi->rows_left) {
		rqi->columns_left = 0;

		int right = _rquad_fixed_to_int(rqi->fx_right_edge);
		rqi->x = _rquad_fixed_to_int(rqi->fx_left_edge);
		if(right <= rqi->x) {
			rqi->columns_left = 0;
		} else {
			rqi->columns_left = right - rqi->x;
		}
		rqi->fs = rqi->s_left - _RQUAD_FIXED_HALF;
		rqi->ft = rqi->t_left - _RQUAD_FIXED_HALF;

		if(rqi->y < 0) {
			// XXX degenerate Opteron specific crap? see mesa original
			rqi->columns_left = 0;
		}

		rqi->rows_left--;
		goto do_col;
	}

	do_subtri:
	if(rqi->sub_triangles_left) {
		_rquad_setup_sub_triangle(rq);
		rqi->sub_triangles_left--;
		goto do_row;
	}

	if(rqi->triangles_left) {
		_rquad_setup_triangle(rq);
		rqi->triangles_left--;
		goto do_subtri;
	}

	return 0;
}

static inline int rquad_x(struct rquad* rq)
{
	return rq->_internal.x;
}

static inline int rquad_y(struct rquad* rq)
{
	return rq->_internal.y;
}

static inline int _rquad_clamp(int i, int min, int max)
{
	if(i < min) return min;
	if(i > max) return max;
	return i;
}

// XXX is clamping really necessary? seems underflow (below 0) is much more
// common than overflow (above max-1). I've not been able to "epsilon" my way
// out of the problem (i.e. by making the source quad slightly smaller than
// requested to eliminate rounding errors near the border). even "big" epsilons
// near half a pixel big don't eliminate the problem, so maybe there are bugs
// in the rasterizer? very big epsilons "solve" it, but breaks rquad in general


static inline int rquad_s(struct rquad* rq)
{
	return _rquad_clamp(_rquad_fixed_to_int(rq->_internal.fs), 0, rq->src_width-1);
}

static inline int rquad_t(struct rquad* rq)
{
	return _rquad_clamp(_rquad_fixed_to_int(rq->_internal.ft), 0, rq->src_height-1);
}

#endif//__RQUAD_H__
