/*
 * ShareLock.h
 *
 *  Created on: 2015年11月4日
 *      Author: guyadong
 */

#ifndef COMMON_SOURCE_CPP_RWLOCK_H_
#define COMMON_SOURCE_CPP_RWLOCK_H_
#include <cstdlib>
#include <cassert>
#include <atomic>
#include <thread>
#include "raii.h"
using namespace std;
namespace gdface {
inline namespace mt{
/*
 * atomic实现读写资源锁,独占写,共享读,禁止复制构造函数和'='赋值操作符
 * WRITE_FIRST为true时为写优先模式,如果有线程等待读取(m_write_wait_count>0)则等待,优先让写线程先获取锁
 * 允许嵌套加锁
 * readLock/Unlock 实现共享的读取加/解锁,线程数不限
 * writeLock/Unlock 实现独占的写入加/解锁,同时只允许一个线程写入，当有线程在读取时，写入线程阻塞，当写入线程执行时，所有的读取线程都被阻塞。
 */
class RWLock {
#define WRITE_LOCK_STATUS -1
#define FREE_STATUS 0
private:
	/* 初始为0的线程id */
	static const  std::thread::id NULL_THEAD;
	const bool WRITE_FIRST;
	/* 用于判断当前是否是写线程 */
	thread::id m_write_thread_id;
	/* 资源锁计数器,-1为写状态，0为自由状态,>0为共享读取状态 */
	atomic_int m_lock_count;
	/* 等待写线程计数器 */
	atomic_uint m_write_wait_count;
public:
	// 禁止复制构造函数
	RWLock(const RWLock&) = delete;
	// 禁止赋值操作符
	RWLock& operator=(const RWLock&) = delete;
	/* 允许移动构造 */
	RWLock(RWLock&& rv):
		WRITE_FIRST(rv.WRITE_FIRST),
		m_write_thread_id(rv.m_write_thread_id),
		m_lock_count(rv.m_lock_count.load()),
		m_write_wait_count(rv.m_write_wait_count.load()){
		// 初始状态的对象才可以移动
		assert(m_write_thread_id==NULL_THEAD&&0==m_lock_count&&0==m_write_wait_count);
	};
	explicit RWLock(bool writeFirst=false)noexcept;//默认为读优先模式
	virtual ~RWLock()=default;
	int readLock()noexcept;
	int readUnlock()noexcept;
	int writeLock()noexcept;
	int writeUnlock()noexcept;
	raii read_guard()const noexcept{
		return gdface::make_raii(*this,&RWLock::readUnlock,&RWLock::readLock);
	}
	raii write_guard()noexcept{
		return gdface::make_raii(*this,&RWLock::writeUnlock,&RWLock::writeLock);
	}
};
}/* namespace mt */
} /* namespace gdface */

#endif /* COMMON_SOURCE_CPP_RWLOCK_H_ */
