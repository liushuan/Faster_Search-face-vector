#ifndef _SSE_NEON_SEARCH__
#define _SSE_NEON_SEARCH__

#ifdef  __GNUC__
#include <arm_neon.h>
#elif _WIN32
#include <immintrin.h>
#define zq_mm256_fmadd_ps _mm256_fmadd_ps
#endif

namespace sse_neon_search {
	/**
		windows dim shoud be [128, 256, 512] , linux dim should be aligned to 8.
	*/
#ifdef  __GNUC__
	inline float _cal_similarity_avx_neon(float* pt1, float* pt2, int dim)
	{
		float sum = 0;
		float32x4_t sum_vec = vdupq_n_f32(0);//, left_vec, right_vec;
		float *A1 = pt1;
		float *B1 = pt2;
		for (int k = 0; k < dim; k += 8)
		{
			A1 += 4;
			B1 += 4;
			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A1), vld1q_f32(B1));
			A1 += 4;
			B1 += 4;
			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A1), vld1q_f32(B1));
		}
		float32x2_t r = vadd_f32(vget_high_f32(sum_vec), vget_low_f32(sum_vec));
		sum += vget_lane_f32(vpadd_f32(r, r), 0);
		return sum;
	}

#elif _WIN32
	inline float _cal_similarity_avx_neon(float* pt1, float* pt2, int dim)
	{
		if (dim == 128) {
			_declspec(align(32)) float q[8];
			__m256 sum_vec = _mm256_mul_ps(_mm256_load_ps(pt1), _mm256_load_ps(pt2));
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 8), _mm256_load_ps(pt2 + 8), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 16), _mm256_load_ps(pt2 + 16), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 24), _mm256_load_ps(pt2 + 24), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 32), _mm256_load_ps(pt2 + 32), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 40), _mm256_load_ps(pt2 + 40), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 48), _mm256_load_ps(pt2 + 48), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 56), _mm256_load_ps(pt2 + 56), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 64), _mm256_load_ps(pt2 + 64), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 72), _mm256_load_ps(pt2 + 72), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 80), _mm256_load_ps(pt2 + 80), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 88), _mm256_load_ps(pt2 + 88), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 96), _mm256_load_ps(pt2 + 96), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 104), _mm256_load_ps(pt2 + 104), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 112), _mm256_load_ps(pt2 + 112), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 120), _mm256_load_ps(pt2 + 120), sum_vec);
			_mm256_store_ps(q, sum_vec);
			float score = q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
			return score;
		}
		else if (dim == 256) {
			_declspec(align(32)) float q[8];
			__m256 sum_vec = _mm256_mul_ps(_mm256_load_ps(pt1), _mm256_load_ps(pt2));
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 8), _mm256_load_ps(pt2 + 8), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 16), _mm256_load_ps(pt2 + 16), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 24), _mm256_load_ps(pt2 + 24), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 32), _mm256_load_ps(pt2 + 32), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 40), _mm256_load_ps(pt2 + 40), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 48), _mm256_load_ps(pt2 + 48), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 56), _mm256_load_ps(pt2 + 56), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 64), _mm256_load_ps(pt2 + 64), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 72), _mm256_load_ps(pt2 + 72), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 80), _mm256_load_ps(pt2 + 80), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 88), _mm256_load_ps(pt2 + 88), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 96), _mm256_load_ps(pt2 + 96), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 104), _mm256_load_ps(pt2 + 104), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 112), _mm256_load_ps(pt2 + 112), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 120), _mm256_load_ps(pt2 + 120), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 128), _mm256_load_ps(pt2 + 128), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 136), _mm256_load_ps(pt2 + 136), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 144), _mm256_load_ps(pt2 + 144), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 152), _mm256_load_ps(pt2 + 152), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 160), _mm256_load_ps(pt2 + 160), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 168), _mm256_load_ps(pt2 + 168), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 176), _mm256_load_ps(pt2 + 176), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 184), _mm256_load_ps(pt2 + 184), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 192), _mm256_load_ps(pt2 + 192), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 200), _mm256_load_ps(pt2 + 200), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 208), _mm256_load_ps(pt2 + 208), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 216), _mm256_load_ps(pt2 + 216), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 224), _mm256_load_ps(pt2 + 224), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 232), _mm256_load_ps(pt2 + 232), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 240), _mm256_load_ps(pt2 + 240), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 248), _mm256_load_ps(pt2 + 248), sum_vec);
			_mm256_store_ps(q, sum_vec);
			float score = q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
			return score;
		}
		else if (dim == 512) {
			_declspec(align(32)) float q[8];
			__m256 sum_vec = _mm256_mul_ps(_mm256_load_ps(pt1), _mm256_load_ps(pt2));
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 8), _mm256_load_ps(pt2 + 8), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 16), _mm256_load_ps(pt2 + 16), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 24), _mm256_load_ps(pt2 + 24), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 32), _mm256_load_ps(pt2 + 32), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 40), _mm256_load_ps(pt2 + 40), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 48), _mm256_load_ps(pt2 + 48), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 56), _mm256_load_ps(pt2 + 56), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 64), _mm256_load_ps(pt2 + 64), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 72), _mm256_load_ps(pt2 + 72), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 80), _mm256_load_ps(pt2 + 80), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 88), _mm256_load_ps(pt2 + 88), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 96), _mm256_load_ps(pt2 + 96), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 104), _mm256_load_ps(pt2 + 104), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 112), _mm256_load_ps(pt2 + 112), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 120), _mm256_load_ps(pt2 + 120), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 128), _mm256_load_ps(pt2 + 128), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 136), _mm256_load_ps(pt2 + 136), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 144), _mm256_load_ps(pt2 + 144), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 152), _mm256_load_ps(pt2 + 152), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 160), _mm256_load_ps(pt2 + 160), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 168), _mm256_load_ps(pt2 + 168), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 176), _mm256_load_ps(pt2 + 176), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 184), _mm256_load_ps(pt2 + 184), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 192), _mm256_load_ps(pt2 + 192), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 200), _mm256_load_ps(pt2 + 200), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 208), _mm256_load_ps(pt2 + 208), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 216), _mm256_load_ps(pt2 + 216), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 224), _mm256_load_ps(pt2 + 224), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 232), _mm256_load_ps(pt2 + 232), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 240), _mm256_load_ps(pt2 + 240), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 248), _mm256_load_ps(pt2 + 248), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 256), _mm256_load_ps(pt2 + 256), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 264), _mm256_load_ps(pt2 + 264), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 272), _mm256_load_ps(pt2 + 272), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 280), _mm256_load_ps(pt2 + 280), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 288), _mm256_load_ps(pt2 + 288), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 296), _mm256_load_ps(pt2 + 296), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 304), _mm256_load_ps(pt2 + 304), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 312), _mm256_load_ps(pt2 + 312), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 320), _mm256_load_ps(pt2 + 320), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 328), _mm256_load_ps(pt2 + 328), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 336), _mm256_load_ps(pt2 + 336), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 344), _mm256_load_ps(pt2 + 344), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 352), _mm256_load_ps(pt2 + 352), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 360), _mm256_load_ps(pt2 + 360), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 368), _mm256_load_ps(pt2 + 368), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 376), _mm256_load_ps(pt2 + 376), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 384), _mm256_load_ps(pt2 + 384), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 392), _mm256_load_ps(pt2 + 392), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 400), _mm256_load_ps(pt2 + 400), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 408), _mm256_load_ps(pt2 + 408), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 416), _mm256_load_ps(pt2 + 416), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 424), _mm256_load_ps(pt2 + 424), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 432), _mm256_load_ps(pt2 + 432), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 440), _mm256_load_ps(pt2 + 440), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 448), _mm256_load_ps(pt2 + 448), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 456), _mm256_load_ps(pt2 + 456), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 464), _mm256_load_ps(pt2 + 464), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 472), _mm256_load_ps(pt2 + 472), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 480), _mm256_load_ps(pt2 + 480), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 488), _mm256_load_ps(pt2 + 488), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 496), _mm256_load_ps(pt2 + 496), sum_vec);
			sum_vec = zq_mm256_fmadd_ps(_mm256_load_ps(pt1 + 504), _mm256_load_ps(pt2 + 504), sum_vec);
			_mm256_store_ps(q, sum_vec);
			float score = q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
			return score;
		}
		else {
			return -1;
		}
	}
#endif
	
}

#endif
