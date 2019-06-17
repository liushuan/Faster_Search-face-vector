/*
 * gf_utility.h
 *
 *  Created on: 2015年11月16日
 *      Author: guyadong
 */

#ifndef COMMON_SOURCE_CPP_GF_UTILITY_H_
#define COMMON_SOURCE_CPP_GF_UTILITY_H_
#include <type_traits>
#include <memory>

#if !defined(_MSC_VER)&&__cplusplus<=201103L
namespace std{
/*
 * 创建普通对象指针unique_ptr
 */
template<typename T,typename ...Args>inline
typename enable_if<!is_array<T>::value,	unique_ptr<T> >::type
make_unique(Args&&...args){
	return unique_ptr<T>(new T(std::forward<Args>(args)...));
}
/*
 * 创建数组对象的unique_ptr
 * T必须是动态数组类型(如: int[])，且不能是定长数组(如 int[20])
 */
template<typename T>inline
typename enable_if<is_array<T>::value&& extent<T>::value==0,unique_ptr<T>>::type
make_unique(size_t size){
	using U=typename remove_extent<T>::type;
	return unique_ptr<T>(new U[size]());
}
template<typename T,typename ...Args>
typename enable_if<extent<T>::value != 0,	void>::type
make_unique(Args&&...)=delete;
}
#endif /* !defined(_MSC_VER)&&__cplusplus<=201103L */
namespace std {
	/* 创建数组类型的shared_ptr */
	template<typename T>
	shared_ptr<T> make_shared_array(size_t size) {
		static_assert(!is_array<T>::value, "T must not be array");
		return shared_ptr<T>(new T[size], default_delete<T[]>());
	}

}
namespace gdface{
/* 修改常量
 * 如:
 * const size_t c=10000;
 * modify_const(c,5ULL);
 */
template <typename T>
void inline modify_const(const T& const_var,const T &new_value)noexcept{
	auto &ref_var =const_cast<T&>(const_var);
	auto &ref_new =const_cast<T&>(new_value);
	ref_var=std::move(ref_new);// 转为右值,以适合比如unique_ptr这种不提供复制操作符的对象
}

/* 判断T有没有==操作符 */
template <typename T>
struct has_equal_operator{
    template<typename U>  static auto test(int)->	decltype(std::declval<U>()==std::declval<U>());
    //template<typename U> static auto test(int)->	decltype(declval<U>().operator==(declval<U>()));
	template<typename U> static void test(...);
    enum{value=std::is_same<decltype(test<T>(0)), bool>::value};
    //通过判断test<T>(0)返回值是否为bool来判断是否有==操作符
};
/* 判断T有没有*操作符 */
template <typename T>
struct has_asterisk_operator{
    template<typename U> static auto test(int)->	decltype(*std::declval<U>());
	template<typename U> static void test(...);
    enum{value=std::is_reference<decltype(test<T>(0))>::value};
    //通过判断test<T>(0)返回值是否为引用来判断是否有*操作符
};
/* 判断T有没有>操作符 */
template <typename T>
struct has_gt_operator{
    template<typename U> static auto test(int)->	decltype(std::declval<U>()>std::declval<U>());
	template<typename U> static void test(...);
	enum{value=std::is_same<decltype(test<T>(0)), bool>::value};
	//通过判断test<T>(0)返回值是否为bool来判断是否有>操作符
};
}/* namespace gdface */
/* 宏函数定义的模板函数，检查T是否有名为's'的成员
 * value 为bool型检查结果
 * type为s成员的类型(value为true是有效)
 */
#define has_member(s) \
template<typename T>\
struct has_member_##s{\
	template <typename _T>static auto check(_T)->typename std::decay<decltype(_T::s)>::type;\
	static void check(...);\
	using type=decltype(check(std::declval<T>()));\
	enum{value=!std::is_void<type>::value};\
};
#endif /* COMMON_SOURCE_CPP_GF_UTILITY_H_ */
