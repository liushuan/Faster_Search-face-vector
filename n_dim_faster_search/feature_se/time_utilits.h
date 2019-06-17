/*
 * time_utility.h
 *
 *  Created on: 2016年4月22日
 *      Author: guyadong
 */

#ifndef COMMON_SOURCE_CPP_TIME_UTILITS_H_
#define COMMON_SOURCE_CPP_TIME_UTILITS_H_
#include <ctime>
namespace gdface{
/*
 * 用于运行时间计算的类
 */
struct time_probe_type{
	clock_t _start;
	/* 开始计时 */
	void start(){
		_start=clock();
	}
	/* 计时结束返回时间(秒) */
	float end(){
		return float(clock()-_start)/CLOCKS_PER_SEC;
	}
};
/* TLS变量用于同一线程内的程序执行时间统计 */
thread_local static time_probe_type tls_clock_probe;
} /* namespace gdface */
#endif /* COMMON_SOURCE_CPP_TIME_UTILITS_H_ */
