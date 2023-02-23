#ifndef WorkPool_H
#define WorkPool_H
#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>
#include <future>

template <class T>
class WorkPool {
  public:
    WorkPool(unsigned int maxThreads);
    WorkPool(unsigned int maxThreads,std::string name);
    ~WorkPool();

    unsigned int initiate();
    bool terminate(bool force = false);
    void wait();

    unsigned int addThreads(unsigned int numThreads);
    void addTask(T f);
    bool addThreadWithTask(T f);

    unsigned int getMaxThreads();
    int numTasks();

  private:
    unsigned int _maxThreads;
    unsigned int _users;

    std::atomic_bool _alive;
    std::atomic_uint _currentThreads;

    std::mutex _tMutex;
    std::vector<std::thread> _threads;

    std::mutex _qMutex;
    std::deque<T> _q;
    std::condition_variable _cv;

    void workLoop();
    std::string _name;
};

#endif /* WorkPool_H */