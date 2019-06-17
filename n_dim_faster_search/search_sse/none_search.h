//#ifndef _NONE_SEARCH__
//#define _NONE_SEARCH__
//
//#include <arm_neon.h>
//
//
////#define vld1q_f32 vld1q_f32
////#define vmlaq_f32 vmlaq_f32
////#define vmlaq_f32 _mm_fmadd_ps
//
//
//
//namespace none_search {
//
//	
//	/*static float dot(const float* A, const float* B, int K)
//	{
//		float sum = 0;
//		float32x4_t sum_vec = vdupq_n_f32(0), left_vec, right_vec;
//		for (int k = 0; k<K; k += 4)
//		{
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A), vld1q_f32(B));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + 4), vld1q_f32(B + 4));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + 8), vld1q_f32(B + 8));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + 16), vld1q_f32(B + 16));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + 3), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//		}
//
//
//		float32x2_t r = vadd_f32(vget_high_f32(sum_vec), vget_low_f32(sum_vec));
//		sum += vget_lane_f32(vpadd_f32(r, r), 0);
//
//		return sum;
//	}*/
//
//	inline   float dot(float* A, float* B, int K)
//	{
//		float sum = 0;
//		float32x4_t sum_vec = vdupq_n_f32(0);//, left_vec, right_vec;
//		for (int k = 0; k<K; k += 8)
//		{
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k), vld1q_f32(B + k));
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A + k + 4), vld1q_f32(B + k + 4));
//			//      sum_vec = vmlap_f32(sum_vec, vld1q_f32(A + k+8), vld1q_f32(B+k+8));
//			//      sum_vec = vmlap_f32(sum_vec, vld1q_f32(A + k+12), vld1q_f32(B+k+12));
//		}
//
//		float32x2_t r = vadd_f32(vget_high_f32(sum_vec), vget_low_f32(sum_vec));
//		sum += vget_lane_f32(vpadd_f32(r, r), 0);
//		return sum;
//
//
//
//
//		float sum = 0;
//		float32x4_t sum_vec = vdupq_n_f32(0);//, left_vec, right_vec;
//		float *A1 = A;
//		float *B1 = B;
//		for (int k = 0; k < K; k += 8)
//		{
//			A1 += 4;
//			B1 += 4;
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A1), vld1q_f32(B1));
//			A1 += 4;
//			B1 += 4;
//			sum_vec = vmlaq_f32(sum_vec, vld1q_f32(A1), vld1q_f32(B1));
//		}
//
//		float32x2_t r = vadd_f32(vget_high_f32(sum_vec), vget_low_f32(sum_vec));
//		sum += vget_lane_f32(vpadd_f32(r, r), 0);
//		return sum;
//
//	}
//
//	/*static float dot(const float* pt1, const float* pt2) {
//		float sum = 0;
//		float32x4_t sum_vec = vdupq_n_f32(0);
//
//	    sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1), vld1q_f32(pt2));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 4), vld1q_f32(pt2 + 4));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 8), vld1q_f32(pt2 + 8));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 12), vld1q_f32(pt2 + 12));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 16), vld1q_f32(pt2 + 16));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 20), vld1q_f32(pt2 + 20));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 24), vld1q_f32(pt2 + 24));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 28), vld1q_f32(pt2 + 28));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 32), vld1q_f32(pt2 + 32));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 36), vld1q_f32(pt2 + 36));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 40), vld1q_f32(pt2 + 40));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 44), vld1q_f32(pt2 + 44));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 48), vld1q_f32(pt2 + 48));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 52), vld1q_f32(pt2 + 52));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 56), vld1q_f32(pt2 + 56));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 60), vld1q_f32(pt2 + 60));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 64), vld1q_f32(pt2 + 64));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 68), vld1q_f32(pt2 + 68));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 72), vld1q_f32(pt2 + 72));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 76), vld1q_f32(pt2 + 76));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 80), vld1q_f32(pt2 + 80));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 84), vld1q_f32(pt2 + 84));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 88), vld1q_f32(pt2 + 88));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 92), vld1q_f32(pt2 + 92));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 96), vld1q_f32(pt2 + 96));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 100), vld1q_f32(pt2 + 100));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 104), vld1q_f32(pt2 + 104));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 108), vld1q_f32(pt2 + 108));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 112), vld1q_f32(pt2 + 112));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 116), vld1q_f32(pt2 + 116));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 120), vld1q_f32(pt2 + 120));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 124), vld1q_f32(pt2 + 124));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 128), vld1q_f32(pt2 + 128));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 132), vld1q_f32(pt2 + 132));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 136), vld1q_f32(pt2 + 136));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 140), vld1q_f32(pt2 + 140));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 144), vld1q_f32(pt2 + 144));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 148), vld1q_f32(pt2 + 148));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 152), vld1q_f32(pt2 + 152));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 156), vld1q_f32(pt2 + 156));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 160), vld1q_f32(pt2 + 160));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 164), vld1q_f32(pt2 + 164));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 168), vld1q_f32(pt2 + 168));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 172), vld1q_f32(pt2 + 172));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 176), vld1q_f32(pt2 + 176));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 180), vld1q_f32(pt2 + 180));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 184), vld1q_f32(pt2 + 184));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 188), vld1q_f32(pt2 + 188));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 192), vld1q_f32(pt2 + 192));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 196), vld1q_f32(pt2 + 196));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 200), vld1q_f32(pt2 + 200));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 204), vld1q_f32(pt2 + 204));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 208), vld1q_f32(pt2 + 208));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 212), vld1q_f32(pt2 + 212));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 216), vld1q_f32(pt2 + 216));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 220), vld1q_f32(pt2 + 220));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 224), vld1q_f32(pt2 + 224));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 228), vld1q_f32(pt2 + 228));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 232), vld1q_f32(pt2 + 232));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 236), vld1q_f32(pt2 + 236));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 240), vld1q_f32(pt2 + 240));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 244), vld1q_f32(pt2 + 244));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 248), vld1q_f32(pt2 + 248));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 252), vld1q_f32(pt2 + 252));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 256), vld1q_f32(pt2 + 256));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 260), vld1q_f32(pt2 + 260));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 264), vld1q_f32(pt2 + 264));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 268), vld1q_f32(pt2 + 268));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 272), vld1q_f32(pt2 + 272));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 276), vld1q_f32(pt2 + 276));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 280), vld1q_f32(pt2 + 280));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 284), vld1q_f32(pt2 + 284));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 288), vld1q_f32(pt2 + 288));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 292), vld1q_f32(pt2 + 292));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 296), vld1q_f32(pt2 + 296));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 300), vld1q_f32(pt2 + 300));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 304), vld1q_f32(pt2 + 304));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 308), vld1q_f32(pt2 + 308));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 312), vld1q_f32(pt2 + 312));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 316), vld1q_f32(pt2 + 316));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 320), vld1q_f32(pt2 + 320));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 324), vld1q_f32(pt2 + 324));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 328), vld1q_f32(pt2 + 328));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 332), vld1q_f32(pt2 + 332));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 336), vld1q_f32(pt2 + 336));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 340), vld1q_f32(pt2 + 340));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 344), vld1q_f32(pt2 + 344));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 348), vld1q_f32(pt2 + 348));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 352), vld1q_f32(pt2 + 352));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 356), vld1q_f32(pt2 + 356));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 360), vld1q_f32(pt2 + 360));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 364), vld1q_f32(pt2 + 364));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 368), vld1q_f32(pt2 + 368));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 372), vld1q_f32(pt2 + 372));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 376), vld1q_f32(pt2 + 376));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 380), vld1q_f32(pt2 + 380));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 384), vld1q_f32(pt2 + 384));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 388), vld1q_f32(pt2 + 388));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 392), vld1q_f32(pt2 + 392));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 396), vld1q_f32(pt2 + 396));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 400), vld1q_f32(pt2 + 400));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 404), vld1q_f32(pt2 + 404));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 408), vld1q_f32(pt2 + 408));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 412), vld1q_f32(pt2 + 412));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 416), vld1q_f32(pt2 + 416));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 420), vld1q_f32(pt2 + 420));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 424), vld1q_f32(pt2 + 424));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 428), vld1q_f32(pt2 + 428));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 432), vld1q_f32(pt2 + 432));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 436), vld1q_f32(pt2 + 436));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 440), vld1q_f32(pt2 + 440));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 444), vld1q_f32(pt2 + 444));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 448), vld1q_f32(pt2 + 448));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 452), vld1q_f32(pt2 + 452));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 456), vld1q_f32(pt2 + 456));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 460), vld1q_f32(pt2 + 460));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 464), vld1q_f32(pt2 + 464));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 468), vld1q_f32(pt2 + 468));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 472), vld1q_f32(pt2 + 472));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 476), vld1q_f32(pt2 + 476));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 480), vld1q_f32(pt2 + 480));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 484), vld1q_f32(pt2 + 484));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 488), vld1q_f32(pt2 + 488));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 492), vld1q_f32(pt2 + 492));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 496), vld1q_f32(pt2 + 496));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 500), vld1q_f32(pt2 + 500));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 504), vld1q_f32(pt2 + 504));
//		sum_vec = vmlaq_f32(sum_vec, vld1q_f32(pt1 + 508), vld1q_f32(pt2 + 508));
//
//		float32x2_t r = vadd_f32(vget_high_f32(sum_vec), vget_low_f32(sum_vec));
//		sum += vget_lane_f32(vpadd_f32(r, r), 0);
//
//		return sum;
//	}*/
//}
//
//
//#endif
