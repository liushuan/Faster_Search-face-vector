/*
* sample_log.h
* 输出带文件名和行号的简单调试日志，使用方式类似java的log4j,用'{}'作占位符
* 调用方式;
*  SAMPLE_LOG("hello,{}",name)
*  Created on: 2016年2月23日
*      Author: guyadong
*/
#ifndef COMMON_SOURCE_CPP_SAMPLE_LOG_H_
#define COMMON_SOURCE_CPP_SAMPLE_LOG_H_
#include <string>
#include <vector>
#include <mutex>
#include <regex>
#include <thread>
#include <iostream>
#include <type_traits>
#include <sstream>
#include "string_utils.h"
#include "file_utilits.h"
namespace gdface {
	namespace log {
		// 模板函数，将value输出到stream
		// 非指针类型参数实现
		template<typename E,
			typename TR = std::char_traits<E>,
			typename T>
			typename std::enable_if<!std::is_pointer<T>::value>::type
			_value_output_stream(std::basic_ostream<E, TR>& stream, const T& value) {
			stream << value;
		}
		// 模板函数，将value输出到stream
		// 指针类型参数实现，value为null时输出字符串‘null’
		template<typename E,
			typename TR = std::char_traits<E>, 
			typename T>
		typename std::enable_if<std::is_pointer<T>::value>::type
			_value_output_stream(std::basic_ostream<E, TR>& stream, const T& value) {
			if (nullptr == value) {
				stream << "null";
			}
			else {
				stream << value;
			}
		}
		// 特化函数
		// 当value为string时转为wstring输出到wostream
		inline void _value_output_stream(std::wostream&stream, const std::string& value) {
			stream << to_wide_string(value);
		}
		// 特化函数
		// 当value为wstring时转为string输出到ostream
		inline void _value_output_stream(std::ostream&stream, const std::wstring& value) {
			stream << to_byte_string(value);
		}
		// 特化函数
		// 当value为wchar_t*时转为string输出到ostream
		inline void _value_output_stream(std::ostream&stream, const wchar_t* value) {
			if (nullptr == value) {
				stream << "null";
			}
			else {
				stream << to_byte_string(value);
			}
		}
		// 特化函数
		// 当value为char*时转为wstring输出到ostream
		inline void _value_output_stream(std::wostream&stream, const char* value) {
			if (nullptr == value) {
				stream << "null";
			}
			else {
				stream << to_wide_string(value);
			}
		}
		/* 终止递归函数 */
		template<typename E,
			typename TR = std::char_traits<E>,
			typename AL = std::allocator<E>>
		void _sm_log_output(std::basic_ostream<E, TR>& stream, const std::vector<std::basic_string<E, TR, AL>>& format, int& idx) {
		}
		// 模板递归处理可变参数，每次处理一个参数
		// T 为第一个参数类型
		template<typename E,
			typename TR = std::char_traits<E>,
			typename AL = std::allocator<E>,
			typename T, typename ...Args>
		void _sm_log_output(std::basic_ostream<E, TR>& stream, const std::vector<std::basic_string<E, TR, AL>>& format, int& idx, const T& first, Args...rest) {
			if (idx < (int)format.size()) {
				_value_output_stream(stream, format[idx]);
				if (idx < (int)format.size() - 1) {
					_value_output_stream(stream, first);
				}
				_sm_log_output(stream, format, ++idx, rest...);
			}
		}
		// 调用递归模板函数_sm_log_output 输出所有可变参数
		// E为基本元素数据类型,支持char,wchar_t，
		// 对应的stream支持ostream,wostream,fromat支持string,wstring
		template<typename E,
			typename TR = std::char_traits<E>,
			typename AL = std::allocator<E>,
			typename _str_type = std::basic_string<E, TR, AL>,
			typename ...Args>
		void sample_log(std::basic_ostream<E, TR>& stream, const char* file, int line, const std::basic_string<E, TR, AL>& format, Args...args) {
			const static char delim[] = "{}";
			static std::once_flag oc;
			std::call_once(oc, [] {
#ifdef _MSC_VER
				std::locale loc(std::locale(), "", LC_CTYPE);
				std::wcout.imbue(loc);
				std::wcerr.imbue(loc);
				std::wclog.imbue(loc);
#elif defined(__GNUC__)
				std::locale::global(std::locale(""));
#endif
			});
			
			// {}为占位符
			auto vf = split(format, std::string("\\{\\}"));
			if (end_with(format, delim)) {
				// 末尾插入空字符串
				vf.push_back(_str_type());
			}
			if (nullptr != file) {
				auto fn = get_file_name(file);
				// 只显示文件名
				auto this_id = std::this_thread::get_id();
				stream << "[" << this_id << "](" << fn.c_str() << ":" << line << ") ";
			}
			int index = 0;
			_sm_log_output(stream, vf, index, args...);
			// 输入参数 少于占位符数目，则原样输出格式化字符串
			for (; index < (int)vf.size(); ++index) {
				stream << vf[index];
				if (index < (int)vf.size() - 1) {
					stream << delim;
				}
			}
			stream << std::endl;
		}
		// 局部特化函数
		// 当format为指针类型时，转为wstring或string
		template<typename E,
			typename TR = std::char_traits<E>,
			typename AL = std::allocator<E>,
			typename ...Args>
			void sample_log(std::basic_ostream<E, TR>& stream, const char* file, int line, const E* format, Args...args) {
			sample_log(stream, file, line, std::basic_string<E, TR, AL>(format), args...);
		}
		// 局部特化函数，
		// 当format为char*类型而stream为wostream类型时，将format转为wstring
		template<typename ...Args>
			void sample_log(std::wostream& stream, const char* file, int line, const char* format, Args...args) {
			sample_log(stream, file, line, to_wide_string(format), args...);
		}
		// 局部特化函数，
		// 当format为wchar_t类型而stream为ostream类型时，将format转为string
		template<typename ...Args>
		void sample_log(std::ostream& stream, const char* file, int line, const wchar_t* format, Args...args) {
			sample_log(stream, file, line, to_byte_string(format), args...);
		}
		// 按 format提供的简单格式输出参数到stream
		template<typename E,
			typename TR = std::char_traits<E>,
			typename AL = std::allocator<E>,
			typename _str_type = std::basic_string<E, TR, AL>,
			typename ...Args>
		void sample_format(std::basic_ostream<E, TR>& stream, const std::basic_string<E, TR, AL>& format, Args...args) {
			sample_log(stream, nullptr, 0, format, args...);
		}
		
		template<typename E,
			typename TR = std::char_traits<E>,
			typename AL = std::allocator<E>,
			typename _str_type = std::basic_string<E, TR, AL>,
			typename ...Args>
			void sample_format(std::basic_ostream<E, TR>& stream, const E* format, Args...args) {
			sample_format(stream,_str_type(format), args...);
		}
		// 按 format提供的简单格式输出参数为string/wstring
		template<typename E,
			typename TR = std::char_traits<E>,
			typename AL = std::allocator<E>,
			typename _str_type = std::basic_string<E, TR, AL>,
			typename ...Args>
		_str_type sample_format(const std::basic_string<E, TR, AL>& format, Args...args) {
			std::basic_stringstream<E,TR> ss;
			sample_format(ss, format, args...);
			return ss.str();
		}
		template<typename E,
			typename TR = std::char_traits<E>,
			typename AL = std::allocator<E>,
			typename _str_type = std::basic_string<E, TR, AL>,
			typename ...Args>
		_str_type sample_format(const E* format, Args...args) {
			std::basic_stringstream<E, TR> ss;
			sample_format(ss, _str_type(format), args...);
			return ss.str();
		}
	} /* namespace log */
// 定义使用 ostream 还是 wostream作为输出流
// 默认使用 wostream 输出,以确保宽字符集信息(如中文)可正确显示
#ifdef _SL_USE_BYTE_STREAM
#	define __SL_STREAM_OUT__ std::cout
#	define __SL_STREAM_ERR__ std::cerr
#	define __SL_STREAM_LOG__ std::clog
#else
#	define __SL_STREAM_OUT__ std::wcout
#	define __SL_STREAM_ERR__ std::wcerr
#	define __SL_STREAM_LOG__ std::wclog
#endif

#define SAMPLE_LOG_STREAM(stream,format,...) ::gdface::log::sample_log(stream,__FILE__,__LINE__,format,##__VA_ARGS__)

// 向std::cout输出带文件名和行号的信息,{}为占位符,调用示例
// SAMPLE_LOG("hello,{} {}","world",2018);
// 输出:hello,world 2018
// NOTE:
// 因为gdface::log::sample_log函数中调用了std::call_once函数,
// 所以在linux下编译时务必要加 -lpthread 选项,否则运行时会抛出异常:
// terminate called after throwing an instance of 'std::system_error'
//		what() : Unknown error - 1
#define SAMPLE_OUT(format,...)	SAMPLE_LOG_STREAM(__SL_STREAM_OUT__,format,##__VA_ARGS__)
#define SAMPLE_ERR(format,...)	SAMPLE_LOG_STREAM(__SL_STREAM_ERR__,format,##__VA_ARGS__)
#define SAMPLE_LOG(format,...)	SAMPLE_LOG_STREAM(__SL_STREAM_LOG__,format,##__VA_ARGS__)

#define SAMPLE_OUT_IF(exp,format,...) if(exp) SAMPLE_OUT(format,__VA_ARGS__);
#define SAMPLE_ERR_IF(exp,format,...) if(exp) SAMPLE_ERR(format,__VA_ARGS__);
#define SAMPLE_LOG_IF(exp,format,...) if(exp) SAMPLE_LOG(format,__VA_ARGS__);

#define RETURN_CHECK_OUT(exp,ret,format,...) if(!(exp)) { SAMPLE_OUT(format,__VA_ARGS__);return ret;}
#define RETURN_CHECK_ERR(exp,ret,format,...) if(!(exp)) { SAMPLE_ERR(format,__VA_ARGS__);return ret;}
#define RETURN_CHECK_LOG(exp,ret,format,...) if(!(exp)) { SAMPLE_LOG(format,__VA_ARGS__);return ret;}

} /* namespace gdface */
#endif /* COMMON_SOURCE_CPP_SAMPLE_LOG_H_ */

