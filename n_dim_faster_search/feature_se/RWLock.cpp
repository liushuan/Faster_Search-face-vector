/*
 * RWLock.cpp
 *
 *  Created on: 2015年11月9日
 *      Author: guyadong
 */
#include "RWLock.h"
namespace gdface {
inline namespace mt{
/* writeFirst==true 写优先模式 */
RWLock::RWLock(bool writeFirst)noexcept:
	WRITE_FIRST(writeFirst),
	m_write_thread_id(),
	m_lock_count(0),
	m_write_wait_count(0){
}
int RWLock::readLock()noexcept {
	// ==时为独占写状态,不需要加锁
	if (this_thread::get_id() != this->m_write_thread_id) {
		int count;
		if (WRITE_FIRST)//写优先模式下,要检测等待写的线程数为0(m_write_wait_count==0)
			do {
				while (m_write_wait_count > 0 || (count = m_lock_count) < 0); //写锁定时等待
			} while (!m_lock_count.compare_exchange_weak(count, count + 1));
		else
			do {
				while ((count = m_lock_count) <0); //写锁定时等待
			} while (!m_lock_count.compare_exchange_weak(count, count + 1));
	}
	return m_lock_count;
}
int RWLock::readUnlock()noexcept {
	// ==时为独占写状态,不需要加锁
	if (this_thread::get_id() != this->m_write_thread_id)
			--m_lock_count;
	return m_lock_count;
}
int RWLock::writeLock()noexcept{
	if (this_thread::get_id() != this->m_write_thread_id){
		++m_write_wait_count;//写等待计数器加1
		// 没有线程读取时(加锁计数器为0)，置为-1加写入锁，否则等待
		for(int zero=FREE_STATUS;!this->m_lock_count.compare_exchange_weak(zero,WRITE_LOCK_STATUS);zero=FREE_STATUS);
		--m_write_wait_count;//获取锁后,计数器减1
		m_write_thread_id=this_thread::get_id();
	}else
		--m_lock_count;// ==时为独占写状态
	return m_lock_count;
}
int RWLock::writeUnlock()noexcept {
	assert(this_thread::get_id() == this->m_write_thread_id&&WRITE_LOCK_STATUS>=m_lock_count);
	if(WRITE_LOCK_STATUS==m_lock_count){
		m_write_thread_id=NULL_THEAD;
		m_lock_count.store(FREE_STATUS);
	}else
		++m_lock_count;
	return m_lock_count;
}
const std::thread::id RWLock::NULL_THEAD;
}/* namespace mt */
} /* namespace gdface */



