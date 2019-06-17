/*
 * memory_ext.h
 *
 *  Created on: 2016年8月21日
 *      Author: guyadong
 */

#ifndef COMMON_SOURCE_CPP_MEMORY_EXT_H_
#define COMMON_SOURCE_CPP_MEMORY_EXT_H_
#include <memory>
namespace std{
template<typename _T>
inline bool operator==(const std::shared_ptr<_T>& __a, const _T &__b) noexcept{ return nullptr==__a?false:*__a== __b; }
template<typename _T>
inline bool operator==(const _T &__a, const std::shared_ptr<_T>& __b) noexcept{ return nullptr==__b?false:*__b== __a; }
template<typename _T>
inline bool operator!=(const std::shared_ptr<_T>& __a, const _T &__b) noexcept{ return nullptr==__a?false:*__a!= __b; }
template<typename _T>
inline bool operator!=(const _T &__a, const std::shared_ptr<_T>& __b) noexcept{ return nullptr==__b?false:*__b!= __a; }
template<typename _T>
inline bool operator>(const std::shared_ptr<_T>& __a, const _T &__b) noexcept{ return nullptr==__a?false:*__a> __b; }
template<typename _T>
inline bool operator>(const _T &__a, const std::shared_ptr<_T>& __b) noexcept{ return nullptr==__b?false:*__b> __a; }
template<typename _T>
inline bool operator>=(const std::shared_ptr<_T>& __a, const _T &__b) noexcept{ return nullptr==__a?false:*__a>= __b; }
template<typename _T>
inline bool operator>=(const _T &__a, const std::shared_ptr<_T>& __b) noexcept{ return nullptr==__b?false:*__b>= __a; }
template<typename _T>
inline bool operator<(const std::shared_ptr<_T>& __a, const _T &__b) noexcept{ return nullptr==__a?false:*__a< __b; }
template<typename _T>
inline bool operator<(const _T &__a, const std::shared_ptr<_T>& __b) noexcept{ return nullptr==__b?false:*__b< __a; }
template<typename _T>
inline bool operator<=(const std::shared_ptr<_T>& __a, const _T &__b) noexcept{ return nullptr==__a?false:*__a<= __b; }
template<typename _T>
inline bool operator<=(const _T &__a, const std::shared_ptr<_T>& __b) noexcept{ return nullptr==__b?false:*__b<= __a; }
} /*namespace std*/
#endif /* COMMON_SOURCE_CPP_MEMORY_EXT_H_ */
