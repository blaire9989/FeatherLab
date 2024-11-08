#pragma once
#include <algorithm>
#include <thread>
#include <functional>
#include <vector>

static void parallel_for(unsigned nb_elements, std::function<void (int start, int end)> functor) {
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);
    nb_threads = nb_threads <= nb_elements ? nb_threads : nb_elements;
    unsigned batch_size = nb_elements / nb_threads;
    unsigned batch_remainder = nb_elements % nb_threads;
    std::vector< std::thread > my_threads(nb_threads);
    for (unsigned i = 0; i < nb_threads; i++) {
        if (i < batch_remainder) {
            int start = i * (batch_size + 1);
            my_threads[i] = std::thread(functor, start, start + batch_size + 1);
        } else {
            int start = batch_remainder * (batch_size + 1) + (i - batch_remainder) * batch_size;
            my_threads[i] = std::thread(functor, start, start + batch_size);
        }
    }
    std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
}