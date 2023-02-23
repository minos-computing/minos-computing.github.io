#include "WorkPool.h"
#include <iostream>

#define DPRINTF(...) fprintf(stderr, __VA_ARGS__)

template <class T>
WorkPool<T>::WorkPool(unsigned int maxThreads) : _maxThreads(maxThreads),
                                                  _users(0),
                                                  _alive(true),
                                                  _currentThreads(0),_name("pool") {
}

template <class T>
WorkPool<T>::WorkPool(unsigned int maxThreads,std::string name) : _maxThreads(maxThreads),
                                                  _users(0),
                                                  _alive(true),
                                                  _currentThreads(0),_name(name) {
}

template <class T>
WorkPool<T>::~WorkPool() {
   std::cout << "[NVDLAEMU] "
              << "deleting thread pool before: "<<_name<<" " << _users << " " << _currentThreads << " " << std::endl;
    terminate(true);
    std::cout << "[NVDLAEMU] "
              << "deleting thread pool: "<<_name<<" " << _users << " " << _currentThreads << " " << std::endl;
}

template <class T>
unsigned int WorkPool<T>::addThreads(unsigned int numThreads) {
    unsigned int threadsToAdd = numThreads;
    std::unique_lock<std::mutex> lock(_tMutex);
    if (_alive.load()) {
        _users++;

        unsigned int currentThreads = _threads.size();
        if (threadsToAdd + currentThreads > _maxThreads)
            threadsToAdd = _maxThreads - currentThreads;

        _currentThreads.fetch_add(threadsToAdd);
        for (unsigned int i = 0; i < threadsToAdd; i++)
            _threads.push_back(std::thread([this] { workLoop(); }));
    }
    lock.unlock();
    return threadsToAdd;
}

template <class T>
unsigned int WorkPool<T>::initiate() {
    return addThreads(_maxThreads);
}

template <class T>
bool WorkPool<T>::terminate(bool force) {
    bool ret = false;
    std::unique_lock<std::mutex> lock(_tMutex);
    if (_users) //So we can join and terminate if there are no users
        _users--;
    if (!_users || force) {
        uint64_t cur_size = _q.size();
        uint64_t timeout_cnt =0;
        while(_q.size() && timeout_cnt < 10){ //let q drain as long as it appears to be making progress otherwise time out
            timeout_cnt++;
            std::this_thread::sleep_for (std::chrono::seconds(1));
            if (cur_size != _q.size()){
                timeout_cnt=0;
                cur_size = _q.size();
            }
        }
        if (timeout_cnt >= 10) {
            std::cerr<<"[NVDLAEMU ERROR] priority thread pool: "<<_name<<" timed out with a non empty queue: "<<_q.size()<<std::endl;
        }
        //This is to deal with the conditional variable
        _alive.store(false);
        while (_currentThreads.load())
            _cv.notify_all();
        //At this point we know the threads have exited
        while (_threads.size()) {
            _threads.back().join();
            _threads.pop_back();
        }
        //if the force is called then we can't reuse the pool
        _alive.store(!force);
        ret = true;
    }
    lock.unlock();
    return ret;
}

template <class T>
void WorkPool<T>::addTask(T f) {
    std::unique_lock<std::mutex> lock(_qMutex);
    _q.push_back(std::move(f));
    lock.unlock();
    _cv.notify_one();
}

template <class T>
void WorkPool<T>::workLoop() {
    T task;
    while (_alive.load()) {
        std::unique_lock<std::mutex> lock(_qMutex);

        //Don't make this a wait or else we will never be able to join
        if (_q.empty()) {
            _cv.wait(lock);
        }

        //Check task since we are waking everyone up when it is time to join
        bool popped = !_q.empty();
        if (popped) {
            task = std::move(_q.front());
            _q.pop_front();
        }

        lock.unlock();

        if (popped) {
            task();
            popped = false;
        }
    }
    //This is the end counter we need to decrement
    _currentThreads.fetch_sub(1);
    if (_q.size() > _currentThreads){
        std::cout<<"[TAZER DEBUG] "<<_name<<" not empty while closing!!!! remaining threads: "<<_currentThreads<<" remaining tasks: "<<_q.size()<<std::endl;
    }
}

template <class T>
void WorkPool<T>::wait() {
    bool full = true;

    while (full) {
        std::unique_lock<std::mutex> lock(_qMutex);
        full = !_q.empty();
        lock.unlock();

        if (full) {
            std::this_thread::yield();
        }
    }
}

template <class T>
unsigned int WorkPool<T>::getMaxThreads() {
    return _maxThreads;
}

template <class T>
bool WorkPool<T>::addThreadWithTask(T f) {
    unsigned int ret = addThreads(1);
    addTask(std::move(f));
    return (ret == 1);
}

template <class T>
int WorkPool<T>::numTasks() {
    return _q.size();
}

template class WorkPool<std::function<void()>>;