#ifndef COMMON_SOURCE_CPP_STRING_UTILS_H_
#define COMMON_SOURCE_CPP_STRING_UTILS_H_
#include <string>
#include <cstdio>
#include <vector>
#include <regex>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <locale>
#include <codecvt>
#ifdef __GNUC__
// 关闭 return std::snprintf(__stream,__n,__format.c_str()); 这行代码产生的警告
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
#endif
namespace gdface {
	/*
	 * 简单模拟java.ang.String.format输出
	 * 除算术类型(int,long,float等)之外只支持std::string输出
	 * */
	class string_utils
	{
		// 普通类型直接返回值
		template<typename T>
		static T cvalue(T&& v){
			return std::forward<T>(v);
		}
		// std::string类型返回c-string指针
		static const char* cvalue( std::string& v){
			return v.c_str();
		}
		// 递归结束
		static int snprintf (char *__stream, size_t __n, const std::string &__format){
			if(__n<=0)return 0;
			// 正则表达式用于获取第一个格式化参数的%控制
			static std::regex fmt_first("^[\\s\\S]*?%[-+ #0]?(\\d+|\\*)?(\\.(?:\\d+|\\*))?(?:hh|h|l|ll|j|z|t|L)?[diuoxXfFeEgGaAcspn]");
			if(std::regex_match(__format,fmt_first)){
				// 实际参数数目与格式化字符串不匹配抛出异常
				throw std::logic_error("invalid format string:missing argument");
			}
			// 调用标准snprintf输出
			return std::snprintf(__stream,__n,__format.c_str());
		}
		template<typename ARG1,typename ...Args>
		static int snprintf (char *__stream, size_t __n, const std::string &__format, ARG1&& first,Args&&...args){
			if(__n<=0)return 0;
			//std::regex fmt_pattern("%[-+ #0]?(\\d+|\\*)?(\\.(?:\\d+|\\*))?(?:hh|h|l|ll|j|z|t|L)?[diuoxXfFeEgGaAcspn]");
			// 正则表达式用于获取第一个格式化参数的%控制
			static std::regex fmt_first("^(?:[^%]|%%)*%[-+ #0]?(\\d+|\\*)?(\\.(?:\\d+|\\*))?(?:hh|h|l|ll|j|z|t|L)?[diuoxXfFeEgGaAcspn]");
			std::smatch m;
			if (!std::regex_search  ( __format, m, fmt_first )){
				// 没有找到%格式控制字符串则抛出异常
				std::cerr<<"extra argument provied to printf\n__format:"<<__format<<std::endl;
				throw std::logic_error("extra argument provied to printf");
			}
			// 调用标准snprintf函数输出第一个参数
			int num=std::snprintf(__stream,__n,m[0].str().c_str(),cvalue(std::forward<ARG1>(first)));
			// 递归调用处理剩余的参数,调整缓冲指针和长度
			return num+snprintf(__stream+num,size_t(__n-num),m.suffix().str(),std::forward<Args>(args)...);
		}
	public:
		template<typename ...Args>
		static std::string format(const std::string &__format,Args&&...args){
			std::vector<char> buffer (size_t(1024));
			int num = snprintf ( buffer.data(), size_t(1024), __format, std::forward<Args>(args)...);
			return std::string(buffer.data(),buffer.data()+num);
		}
	};
	inline std::wstring to_wide_string(const std::string& input)
	{
		std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
		return converter.from_bytes(input);
	}
	inline std::string to_byte_string(const std::wstring& input)
	{
		//std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
		std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
		return converter.to_bytes(input);
	}
	inline std::string trim(const std::string &s) {
		if (s.empty()) {
			return const_cast<std::string&>(s);
		}
		std::string out(s);
		out.erase(0, out.find_first_not_of(" "));
		out.erase(out.find_last_not_of(" ") + 1);
		return std::move(out);
	}
	inline bool has_space(const std::string &s) {
		return s.find_first_of(' ') != std::string::npos;
	}
	template<typename E,
		typename TR = std::char_traits<E>,
		typename AL = std::allocator<E>>
	inline std::basic_string<E, TR, AL> toupper(const std::basic_string<E, TR, AL>&src) {
		auto dst = src;
		static const std::locale loc("");
		transform(src.begin(), src.end(), dst.begin(), [&](E c)->E {return std::toupper(c, loc); });
		return dst;
	}
	template<typename E,
		typename TR = std::char_traits<E>,
		typename AL = std::allocator<E>>
	inline std::basic_string<E, TR, AL> tolower(const std::basic_string<E, TR, AL>&src) {
		auto dst = src;
		static const std::locale loc("");
		transform(src.begin(), src.end(), dst.begin(), [&](E c)->E {return std::tolower(c, loc); });
		return dst;
	}
	template<typename E,
		typename TR = std::char_traits<E>,
		typename AL = std::allocator<E>>
	inline bool end_with(const std::basic_string<E, TR, AL>&src, const std::basic_string<E, TR, AL> &suffix) {
		if (src.size() < suffix.size()) {
			return false;
		}
		return src.substr(src.size() - suffix.size()) == suffix;
	}
	inline bool end_with(const std::wstring&src, const std::string&suffix) {
		return end_with(src, to_wide_string(suffix));
	}
	template<typename E,
		typename TR = std::char_traits<E>,
		typename AL = std::allocator<E>>
	inline bool end_with(const std::basic_string<E, TR, AL>&src, const E*suffix) {
		return end_with(src, std::basic_string<E, TR, AL>(suffix));
	}
	template<typename E,
		typename TR = std::char_traits<E>,
		typename AL = std::allocator<E>>
	inline bool start_with(const std::basic_string<E, TR, AL>&src, const std::basic_string<E, TR, AL> &prefix) {
		return src.substr(0, prefix.size()) == prefix;
	}
	inline bool start_with(const std::wstring&src, const std::string&suffix) {
		return start_with(src, to_wide_string(suffix));
	}
	template<typename E,
		typename TR = std::char_traits<E>,
		typename AL = std::allocator<E>>
		inline bool start_with(const std::basic_string<E, TR, AL>&src, const E*suffix) {
		return start_with(src, std::basic_string<E, TR, AL>(suffix));
	}
	template<typename E,
		typename TR = std::char_traits<E>,
		typename AL = std::allocator<E>,
		typename _str_type = std::basic_string<E, TR, AL>>
	std::vector<_str_type> split(const std::basic_string<E, TR, AL>& in, const std::basic_string<E, TR, AL>& delim) {
		std::basic_regex<E> re{ delim };
		return std::vector<_str_type> {
			std::regex_token_iterator<typename _str_type::const_iterator>(in.begin(), in.end(), re, -1),
				std::regex_token_iterator<typename _str_type::const_iterator>()
		};
	}
	inline std::vector<std::wstring> split(const std::wstring& in, const std::string& delim) {
		return split(in, to_wide_string(delim));
	}
	template<typename E,
		typename TR = std::char_traits<E>,
		typename AL = std::allocator<E>,
		typename _str_type = std::basic_string<E, TR, AL>>
		std::vector<_str_type> split(const std::basic_string<E, TR, AL>& in, const E* delim) {
		return split(in, _str_type(delim));
	}
	// c string版本
	inline std::vector<std::string> split(const char* in, const char* delim) {
		std::regex re{ delim };
		return std::vector<std::string> {
			std::cregex_token_iterator(in, in + strlen(in),re, -1),
			std::cregex_token_iterator()
		};
	}
	// 支持wchar_t宽字符集的版本
	inline std::vector<std::wstring> split(const wchar_t* in, const wchar_t* delim) {
		std::wregex re{ delim };
		return std::vector<std::wstring> {
			std::wcregex_token_iterator(in, in + wcslen(in),re, -1),
			std::wcregex_token_iterator()
		};
	}
} /* namespace gdface */
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#endif /* COMMON_SOURCE_CPP_STRING_UTILS_H_ */
