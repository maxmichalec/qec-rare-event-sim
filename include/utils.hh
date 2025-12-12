#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <random>
#include <chrono>
#include <iostream>

#define DEBUG 0

using namespace std;

// --- Stats ---
struct Stats {
    chrono::milliseconds total;
    chrono::milliseconds mc_seed;
    chrono::milliseconds stim_py_time;
    vector<chrono::milliseconds> metropolis_splits;
    vector<int> stim_sim_calls;
};

inline void summarize_stats(Stats stats) {
    cout << "--- Statistics Summary ---" << endl;
    cout << "Total runtime: " << stats.total.count() / 1000.0 << "s" << endl;
    cout << "Initial MC seed runtime: " << stats.mc_seed.count() / 1000.0 << "s" << endl;
    cout << "Total Stim & PyMatching runtime during Metropolis: " << stats.stim_py_time.count() / 1000.0 << "s" << endl;
    int64_t metropolis_ms = 0.0;
    for (auto split : stats.metropolis_splits) {
        metropolis_ms += split.count();
    }
    cout << "Total Metropolis iters runtime: " << metropolis_ms / 1000.0 << "s" << endl;
    cout << "Average Metropolis iter runtime: " << (metropolis_ms / 1000.0) / stats.metropolis_splits.size() << "s" << endl;
    int stim_calls = 0;
    for (int split : stats.stim_sim_calls) {
        stim_calls += split;
    }
    cout << "Total calls to Stim simulation: " << stim_calls << endl;
    cout << "--- End Summary ---\n" << endl;
};

// --- RNG ---
mt19937 rng(random_device{}());

struct RNG {
    static double uniform() {
        return uniform_real_distribution<double>(0.0, 1.0)(rng);
    }

    static bool bernoulli(double p) {
        return uniform() < p;
    }

    template<typename T>
    static T random_choice(const vector<T>& vec) {
        uniform_int_distribution<int> dist(0, vec.size() - 1);
        return vec[dist(rng)];
    }
};

#endif /* UTILS_H_ */
