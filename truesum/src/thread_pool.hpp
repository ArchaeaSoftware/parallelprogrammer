// A fixed set of workers, parked between calls.
//
// Both accumulation matrices want the same thing: run one function over a
// partition of the columns and wait. Creating threads per call is not viable --
// a batch is often only a hundred microseconds, and thread creation is a real
// fraction of that, measured at 2.73 Gelem/s against 4.64 for a pool.
//
// The calling thread takes slot 0 and runs its own share, so `n` threads means
// n-1 workers and no handoff for the single-threaded case.
#pragma once

#include <condition_variable>
#include <cstddef>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace truesum {
namespace detail {

class ThreadPool {
public:
    explicit ThreadPool(unsigned n) : count_(n < 1 ? 1 : n)
    {
        for (unsigned w = 1; w < count_; ++w) {
            workers_.emplace_back([this, w] {
                unsigned seen = 0;
                for (;;) {
                    std::unique_lock<std::mutex> lk(m_);
                    go_.wait(lk, [&] { return stop_ || epoch_ != seen; });
                    if (stop_) return;
                    seen = epoch_;
                    const std::function<void(unsigned)> fn = job_;
                    lk.unlock();
                    fn(w);
                    lk.lock();
                    if (0 == --outstanding_) done_.notify_one();
                }
            });
        }
    }

    ~ThreadPool()
    {
        {
            std::lock_guard<std::mutex> lk(m_);
            stop_ = true;
        }
        go_.notify_all();
        for (auto &t : workers_) t.join();
    }

    ThreadPool(const ThreadPool &) = delete;
    ThreadPool &operator=(const ThreadPool &) = delete;

    unsigned count() const { return count_; }

    // Runs fn(0) on this thread and fn(w) on each worker, returning when all
    // have finished.
    void run(const std::function<void(unsigned)> &fn)
    {
        if (workers_.empty()) {
            fn(0);
            return;
        }
        {
            std::lock_guard<std::mutex> lk(m_);
            job_ = fn;
            outstanding_ = static_cast<unsigned>(workers_.size());
            ++epoch_;
        }
        go_.notify_all();
        fn(0);
        std::unique_lock<std::mutex> lk(m_);
        done_.wait(lk, [&] { return 0 == outstanding_; });
    }

    // The [begin, end) of `total` items belonging to slot `slot`, in
    // contiguous blocks so a worker's items stay near one another in memory.
    void partition(std::size_t &begin, std::size_t &end, std::size_t total,
                   unsigned slot) const
    {
        const std::size_t per = (total + count_ - 1) / count_;
        begin = per * slot;
        if (begin > total) begin = total;
        end = begin + per;
        if (end > total) end = total;
    }

private:
    unsigned count_;
    std::vector<std::thread> workers_;
    std::mutex m_;
    std::condition_variable go_, done_;
    std::function<void(unsigned)> job_;
    unsigned epoch_ = 0;
    unsigned outstanding_ = 0;
    bool stop_ = false;
};

}  // namespace detail
}  // namespace truesum
