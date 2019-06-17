/*
 * common_utilits.h
 *
 *  Created on: 2016年1月7日
 *      Author: guyadong
 */

#ifndef COMMON_SOURCE_CPP_COMMON_UTILITS_H_
#define COMMON_SOURCE_CPP_COMMON_UTILITS_H_
#include <cstring>
#include <string>
#include <memory>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <stdexcept>
using namespace std;
namespace gdface {
inline  namespace com_utilits{
inline char hex_to_byte(char c){
	if(c>='0' && c<='9')return c-'0';
	else if(c>='a' && c<='f')return c-'a'+10;
	else if(c>='A' && c<='F')return c-'A'+10;
	else
		throw invalid_argument("invalid hex char[0-9a-fA-F]");
}
inline void hex_to_bytes(const char* hex_str, void *out, size_t size) {
	static const string err_invalid_hex_str="the argument hex_str is not a valid hex string";
	if (nullptr==hex_str||nullptr == out)
		throw std::invalid_argument("the argument 'hex_str' and 'bytes' must not be null");
	auto len=strlen(hex_str);
	if (0==len||(len&1))
		throw std::invalid_argument(err_invalid_hex_str);
	auto limit = std::min(len >> 1, size) << 1;
	try{
		auto bytes=reinterpret_cast<uint8_t*>(out);
		for (auto i = 0; i < limit; i += 2) {
			bytes[i >> 1] = (hex_to_byte(hex_str[i]) << 4) + hex_to_byte(hex_str[i + 1]);
		}
	}catch(invalid_argument &e){
		throw std::invalid_argument(string(err_invalid_hex_str).append(" because:").append(e.what()));
	}
}

inline vector<uint8_t> hex_to_bytes(const char* hex_str) {
	vector<uint8_t> v(strlen(hex_str) >> 1);
	hex_to_bytes(hex_str, v.data(), v.size());
	return std::move(v);
}
inline vector<uint8_t> hex_to_bytes(const string &hex_str) {
	return hex_to_bytes(hex_str.data());
}
inline void bytes_to_hex(const void *in, size_t size, char* hex_str, size_t limit) {
	const static char hex_char[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
	if (nullptr == in || nullptr == hex_str)
		throw std::invalid_argument("the argument 'bytes' and 'hex_str' must not be null");
	if (limit < size << 1)
		throw std::invalid_argument("the capacity of out buffer  'hex_str' can't fit size of input ");
	auto bytes=reinterpret_cast<const uint8_t*>(in);
	for (decltype(size) i = 0; i < size; ++i) {
		hex_str[i << 1] = hex_char[bytes[i] >> 4]; // 高四位
		hex_str[(i << 1) + 1] = hex_char[bytes[i] & 0x0f]; // 低四位
	}
}
inline vector<char> bytes_to_hex(const void *bytes, size_t size) {
	vector<char> v((size << 1) + 1);
	bytes_to_hex(bytes, size, v.data(), v.size());
	v[size << 1] = 0;
	return std::move(v);
}
inline vector<char> bytes_to_hex(const vector<uint8_t> &v) {
	return bytes_to_hex(v.data(), v.size());
}
inline string bytes_to_hex_string(const void *bytes, size_t size) {
	return string(bytes_to_hex(bytes, size).data());
}
inline string bytes_to_hex_string(const vector<uint8_t> &v) {
	return bytes_to_hex_string(v.data(), v.size());
}

} /* namespace com_utilits */

}  /* namespace gdface */

#endif /* COMMON_SOURCE_CPP_COMMON_UTILITS_H_ */
