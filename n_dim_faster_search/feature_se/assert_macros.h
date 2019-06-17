/*
 * assert_utilits.h
 *
 *  Created on: 2016年2月3日
 *      Author: 10km
 */

#ifndef COMMON_SOURCE_CPP_ASSERT_MACROS_H_
#define COMMON_SOURCE_CPP_ASSERT_MACROS_H_
#include <stdexcept>
#include <string>
#define _DEF_STRING(x) #x
#define DEF_TO_STRING(x) _DEF_STRING(x)
#define SOURCE_AT __FILE__ ":" DEF_TO_STRING(__LINE__)
#define ERROR_STR(msg) std::string(SOURCE_AT ":").append(msg)

#define throw_except_if_msg(except,expression,msg) \
	if(expression)\
		throw except(ERROR_STR(msg));
#define throw_except_if(except,expression) throw_except_if_msg(except,expression,#expression)
#define throw_if_msg(expression,msg) throw_except_if_msg(std::invalid_argument,expression,msg)
#define throw_if(expression) throw_except_if(std::invalid_argument,expression)
#define throw_except_if_null(except,p) throw_except_if_msg(except,nullptr==p,#p" is null")
#define throw_if_null(p) throw_if_msg(nullptr==p,#p" is null")

#endif /* COMMON_SOURCE_CPP_ASSERT_MACROS_H_ */
