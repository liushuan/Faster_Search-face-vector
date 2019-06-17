/*
 * file_utilits.h
 *
 *  Created on: 2016年2月23日
 *      Author: guyadong
 */

#ifndef COMMON_SOURCE_CPP_FILE_UTILITS_H_
#define COMMON_SOURCE_CPP_FILE_UTILITS_H_
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <dirent.h>
#ifdef _WIN32
#include<direct.h>
#else
#include <unistd.h>
#endif
#include "raii.h"
#include "assert_macros.h"
namespace gdface{
/* 将二进制写入指定的文件，filename为nullptr或写入失败则抛出异常invalid_argument  */
inline bool save_binary_file(const std::string&filename, const void *bin_data, size_t size) {
	throw_if_null(bin_data)
	gdface::raii_var<std::ofstream> raii_ofs(
			[&]() {return std::ofstream(filename, std::ofstream::binary); },
			[](std::ofstream&ofs) {ofs.close(); });
	throw_if_msg(!*raii_ofs, std::string("can not open file:").append(filename))
		
	raii_ofs->write((const char*)bin_data, size);
#ifndef NDEBUG
	std::cout << filename << " saved,size=" << size << std::endl;
#endif
	return true;
}
inline bool save_binary_file(const char *filename,const void *bin_data, size_t size){
	throw_if_null(filename)
	return save_binary_file(std::string(filename),bin_data,size);
}
/* 以二进制读取指定的istream */
inline std::vector<uint8_t> load_binary_stream(std::istream& is) {
	std::vector<uint8_t> bin_data;
	// get length of file:
	is.seekg(0, is.end);
	auto length = is.tellg();
	is.seekg(0, is.beg);
	bin_data = std::vector<uint8_t>(size_t(length));
	// read data as a block:
	is.read((char*)bin_data.data(), bin_data.size());
	return std::move(bin_data);
}
/* 以二进制读取指定的文件 */
inline std::vector<uint8_t> load_binary_file(const std::string&filename) {
	gdface::raii_var<std::ifstream> raii_ifs(
		[&]() {return std::ifstream(filename, std::ifstream::binary); },
		[](std::ifstream&ifs) {ifs.close(); });
	throw_if_msg(!*raii_ifs, std::string("can not open file:").append(filename))
		return load_binary_stream(*raii_ifs);
}
/* 以二进制读取指定的文件，filename为nullptr或打开失败则抛出异常invalid_argument  */
inline std::vector<uint8_t> load_binary_file(const char *filename){
	throw_if_msg(nullptr==filename,"filename is nullptr")
	return load_binary_file(std::string(filename));
}

/* 以指定的文本文件转为std::string */
inline std::string load_string(const std::string&filename) {
	auto bin = gdface::load_binary_file(filename);
	return std::string((char*)bin.data(), bin.size());
}
/* 以指定的文本文件转为std::string，filename为nullptr或打开失败则抛出异常invalid_argument  */
inline std::string load_string(const char *filename){
	return load_string(std::string(filename));
}
/* 获取当前路径 */
inline std::string getcwd(){
	raii_var<char*> raii_str(
#ifdef _MSC_VER
			[](){return _getcwd(nullptr,0);},
#else
			[]() {return ::getcwd(nullptr, 0); },
#endif
			[](char*p){free(p);}
			);
	return std::string(*raii_str);
}

#ifdef _WIN32
inline char file_sepator(){
	return '\\';
}
inline const std::string& valid_file_sepator() {
	static const std::string se("\\/");
	return se;
}
#else
inline char file_sepator(){
	return '/';
}
inline const std::string& valid_file_sepator() {
	static const std::string se("/");
	return se;
}
#endif

inline bool end_with_file_sepator(const std::string &path) {
	return !path.empty() && std::string::npos != valid_file_sepator().find(path[path.size() - 1]);
}
inline bool start_with_file_sepator(const std::string &path) {
	return !path.empty() && std::string::npos != valid_file_sepator().find(path[0]);
}
inline std::string path_concate(const std::string &path1,const std::string &path2){
	return path1.empty()||end_with_file_sepator(path1) 
		? path1 + path2 
		: path1 + file_sepator() + path2;
}

inline bool is_folder(const char* dir_name){
	throw_if(nullptr==dir_name);
	auto dir =opendir(dir_name);
	if(dir){
		closedir(dir);
		return true;
	}
	return false;
}
inline bool is_folder(const std::string &dir_name){
	throw_if(dir_name.empty());
	return is_folder(dir_name.data());
}
using file_filter_type=std::function<bool(const char*,const char*)>;
const static  file_filter_type default_ls_filter=[](const char*,const char*){return true;};
/*
 * 列出指定目录的所有文件(不包含目录)执行，对每个文件执行filter过滤器
 * sub为true时为目录递归
 * 返回每个文件的全路径名
*/
inline  std::vector<std::string> for_each_file(const std::string&dir_name,file_filter_type filter,bool sub=false){
	std::vector<std::string> files;
	auto dir =opendir(dir_name.data());
	struct dirent *ent;
	if(dir){
		raii guard([&](){closedir(dir);});
		while ((ent = readdir (dir)) != nullptr) {
			auto p = std::string(dir_name).append({ file_sepator() }).append(ent->d_name);
			if(sub){
				if ( 0== strcmp (ent->d_name, "..") || 0 == strcmp (ent->d_name, ".")){
					continue;
				}else if(is_folder(p)){
					auto r= for_each_file(p,filter,sub);
					files.insert(files.end(),r.begin(),r.end());
					continue;
				}
			}
			if (sub||!is_folder(p))
				if(filter(dir_name.data(),ent->d_name))
					files.emplace_back(p);
		}
	}
	return files;
}
/*
 * 列出指定目录的所有文件
 * sub为true时为目录递归
 * 返回每个文件的全路径名
 */
inline std::vector<std::string> ls(const std::string&dir_name, bool sub = false) {
	return for_each_file(dir_name, default_ls_filter, sub);
}

/*
 * 从路径中获取文件名
*/
inline std::string get_file_name(const std::string &path){
	auto find=path.find_last_of(valid_file_sepator());
	return std::string::npos==find?path:path.substr(find+1);
}
/*
 * 从文件名中查找最后一个.符号，获取文件后缀，
 * 如果找不到或.位于字符串起始位置或位于字符结尾，则返回空字符串
*/
inline std::string get_file_suffix(const std::string &file_name){
	auto find=file_name.find_last_of('.');
	return std::string::npos==find||0==find||file_name.size()-1==find?"":file_name.substr(find+1);
}

/* 终止递归函数 */
inline void args_print(std::ostream& steam){}
/* 使用可变参数模板实现参数打印到输出流(ostream)  */
template<typename T,typename ...Args>
inline void args_print(std::ostream& steam,T first,Args...rest){
	steam<< first;
	args_print(steam,rest...);
}
} /* namespace gdface*/



#endif /* COMMON_SOURCE_CPP_FILE_UTILITS_H_ */
