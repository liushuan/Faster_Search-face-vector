/*
 * intrinsic_wrapper.h
 *
 *  Created on: 2015年10月5日
 *      Author: guyadong
 */


#ifndef COMMON_SOURCE_CPP_INTRINSIC_WRAPPER_H_
#define COMMON_SOURCE_CPP_INTRINSIC_WRAPPER_H_
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
#if _MSC_VER
#include <intrin.h>
/* msc 编译时使用内建函数 (intrinsics) */
#if _M_AMD64 || _M_ARM
int8_t _inline _bsr_int64_(uint64_t num) {
	unsigned long index;
	return _BitScanReverse64(&index, num)?(int8_t)index:-1;
}
#define _BSR_INT64_DEFINED
#endif
int8_t _inline _bsr_int32_(uint32_t num) {
	unsigned long index;
	return _BitScanReverse(&index, num)?(int8_t)index:-1;
}
#elif (__x86_64__||_X86_)  && __GNUC__
#ifdef __x86_64__
inline int8_t _bsr_int64_(uint64_t num) {
	uint64_t count;
	__asm__(
			"bsrq %1, %0\n\t"
			"jnz 1f\n\t"
			"movq $-1,%0\n\t"
			"1:"
			:"=q"(count):"q"(num));
	return count;
}
#define _BSR_INT64_DEFINED
#endif
inline int8_t _bsr_int32_(uint32_t num) {
	uint32_t count;
	__asm__(
			"bsrl %1, %0\n\t"
			"jnz 1f\n\t"
			"movl $-1,%0\n\t"
			"1:"
			:"=r"(count):"r"(num));
	return count;
}
#elif __GNUC__
//gcc 编译时使用内建函数(builtin)
inline int8_t _bsr_int64_(uint64_t num) {
	return num==0?-1:(sizeof(num)<<3)-1-__builtin_clzll(num);
}
inline int8_t _bsr_int32_(uint32_t num) {
	return num==0?-1:(sizeof(num)<<3)-1-__builtin_clz(num);
}
#define _BSR_INT64_DEFINED
#else
_inline int8_t _bsr_int64_(uint64_t num) {
	__int8 count=(sizeof(num)<<3)-1;
	for(uint64_t mask=1LLU<<count;!(num&mask)&&count>=0;count--,mask>>=1);
	return count;
}
_inline int8_t _bsr_int32_(uint32_t num) {
	return _bsr_int64_(num);
}
#define _BSR_INT64_DEFINED
#endif

#ifndef _BSR_INT64_DEFINED
int8_t _inline _bsr_int64_(uint64_t num) {
	int8_t r=_bsr_int32_((uint32_t)(num>>32));
	return r>=0?r+32:_bsr_int32_((uint32_t)(num&0xffffffff));
}
#else
#undef _BSR_INT64_DEFINED
#endif
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
#if _MSC_VER && (_M_AMD64||_M_IX86)
#define __builtin_popcount __popcnt
#define __POPCOUNT_ENABLE	1
#if _M_AMD64
#define __builtin_popcountll __popcnt64
#define __POPCOUNT64_ENABLE	1
#endif
#elif __GNUC__
#define __POPCOUNT_ENABLE	1
#define __POPCOUNT64_ENABLE	1
#endif
namespace gdface {
namespace intrinsics {
	inline int8_t _bsr(uint64_t num) {
		return _bsr_int64_(num);
	}
	inline int8_t _bsr(uint32_t num) {
		return _bsr_int32_(num);
	}
	inline int8_t _bsr(uint16_t num) {
		return _bsr_int32_(num);
	}
	inline int8_t _bsr(uint8_t num) {
		return _bsr_int32_(num);
	}
#if __POPCOUNT64_ENABLE
	inline int8_t _popcount(uint64_t num){
		return (int8_t)__builtin_popcountll(num);
	}
#endif
#if __POPCOUNT_ENABLE
	inline int8_t _popcount(uint32_t num){
		return (int8_t)__builtin_popcount(num);
	}
	inline int8_t _popcount(uint16_t num){
		return (int8_t)__builtin_popcount(num);
	}
	inline int8_t _popcount(uint8_t num){
		return (int8_t)__builtin_popcount(num);
	}
#endif

}/*namespace intrinsics*/
} /* namespace gdface */
#ifdef __builtin_popcount
#undef __builtin_popcount
#endif
#ifdef __builtin_popcountll
#undef __builtin_popcountll
#endif
#ifdef __POPCOUNT_ENABLE
#undef __POPCOUNT_ENABLE
#endif
#ifdef __POPCOUNT64_ENABLE
#undef __POPCOUNT64_ENABLE
#endif
#endif
#endif /* COMMON_SOURCE_CPP_INTRINSIC_WRAPPER_H_ */


