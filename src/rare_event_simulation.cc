#include <vector>
#include <random>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <iterator>
#include <chrono>
#include <nlohmann/json.hpp>

#include "circuit.hh"
#include "utils.hh"

using namespace std;
using namespace nlohmann;

// Initialize global utilities
Stats stats;

// --- Events ---
typedef vector<Fault> Event;

void write_faults_file(const Circuit& circt, const Event& E, const string& fpath) {
    json j;
    j["d"] = circt.d;

    // int n_ancilla = (circt.d * circt.d) - 1;
    for (int r = 0; r < circt.d; ++r) {
        map<int, string> prep_faults;
        vector<map<int, string>> cnot_faults(4);
        map<int, string> meas_faults;
        for (int i = 0; i < circt.gates_per_round; ++i) {
            int g = r * circt.gates_per_round + i;
            Gate gate = circt.gates[g];
            if ((gate.g_type == GateType::INIT_X) || (gate.g_type == GateType::INIT_Z)) {
                if (E[g] != Fault::NONE) {
                    prep_faults.insert(pair<int, string>(point_to_idx(gate.qubit, circt.d), faultToString.at(E[g])));
                }
            } else if (gate.g_type == GateType::CNOT) {
                if (E[g] == Fault::NONE) {
                    // cnot_faults[gate.target].insert(pair<int, string>(point_to_idx(gate.qubit, circt.d), "II"));
                } else {
                    cnot_faults[gate.target].insert(pair<int, string>(point_to_idx(gate.qubit, circt.d), faultToString.at(E[g])));
                }
            } else if ((gate.g_type == GateType::MEASURE_X) || (gate.g_type == GateType::MEASURE_Z)) {
                if (E[g] != Fault::NONE) {
                    meas_faults.insert(pair<int, string>(point_to_idx(gate.qubit, circt.d), faultToString.at(E[g])));
                }
            }
        }
        string round = to_string(r);
        for (auto pair : prep_faults) {
            j["faults"][round]["PREP"][to_string(pair.first)] = pair.second;
        }
        j["faults"][round]["CNOT"] = json::array();
        for (int i = 0; i < 4; ++i) {
            json cnot_seq_json;
            for (auto pair : cnot_faults[i]) {
                cnot_seq_json[to_string(pair.first)] = pair.second;
            }
            j["faults"][round]["CNOT"].push_back(cnot_seq_json);
        }
        for (auto pair : meas_faults) {
            j["faults"][round]["MEASURE"][to_string(pair.first)] = pair.second;
        }
    }

    ofstream of(fpath);
    if (!of.is_open()) {
        exit(1);
    }
    of << j.dump(2) << endl;
    of.close();
}

// --- Decoder ---
bool decode(const Circuit& circt, const Event& E) {
    string faults_fpath = "/tmp/surface_code_faults.json";

    write_faults_file(circt, E, faults_fpath);

    ostringstream oss;
    oss << "./.venv/bin/python3 src/python/sim_stim_faults.py " << faults_fpath;
    string cmd = oss.str();

    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        exit(1);
    }

    char buffer[128];
    string result = "";
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }

    pclose(pipe);

    int causes_logical_error = stoi(result);

#if DEBUG
    if (causes_logical_error != 0) {
        cout << "Found error!" << endl;
    }
#endif

    return (causes_logical_error != 0);
}

// --- Find Seed Logical Error Event ---
// Utility function to find an Event that causes a logical failure to seed Metropolis routine
Event find_logical_error(const Circuit& circt) {
    Event E = Event(circt.G, Fault::NONE);

    while (true) {
        for (int g = 0; g < circt.G; ++g) {
            if (RNG::bernoulli(0.1)) {
                int f_idx = uniform_int_distribution<int>(0, circt.gates[g].p_faults.size() - 1)(rng);
                auto it = circt.gates[g].p_faults.begin();
                advance(it, f_idx);
                Fault f = it->first;
                E[g] = f;
            }
        }

        if (decode(circt, E)) {
            cout << "Found seed error." << endl;
            return E;
        }

        // Reset event list
        for (int i = 0; i < circt.G; ++i) {
            E[i] = Fault::NONE;
        }
    }
}

// --- Rare Event Simulation Components ---
double log_event_probability(const Circuit& circt, const Event& E, double p) {
    double log_p_E = 1.0;
    for (int i = 0; i < circt.G; ++i) {
        double p_g = (circt.gates[i].p_failure > 0.0) ? circt.gates[i].p_failure : p;  // probability that gate failed
        Fault f = E[i];
        if (f != Fault::NONE) {
            double p_failure = circt.gates[i].p_faults.at(f);
            log_p_E += log(p_g) + log(p_failure);
        } else {  // gate didn't have fault
            log_p_E += log(1.0 - p_g);
        }
    }
    return log_p_E;
}

inline double log_event_probability_ratio(const Circuit& circt, const Event& E, double p_i, double p_j) {
    return log_event_probability(circt, E, p_i) - log_event_probability(circt, E, p_j);
}

// - Splitting Sequence -
double splitting_sequence_next(const Circuit& circt, double p_i) {
    double w_i = max((double)(circt.d / 2.0), p_i * circt.G);
    double p_next = p_i * pow(2.0, -1.0 / sqrt(w_i));
    return p_next;
}

// double splitting_sequence_next(const Circuit& circt, double p_i) {
//     // More fine grained sequence
//     return p_i * 0.9;
// }

// - Metropolis Routine -
vector<Event> run_metropolis_routine(const Circuit& circt, double p, const Event& E_init, int N) {
    auto start = chrono::high_resolution_clock::now();

    vector<Event> failing_events;
    Event E = E_init;

    // Sample (gate, fault) pairs until N failing events
    int iters = 0;
    // int stim_sim_calls = 0;
    while ((int) failing_events.size() < N) {
        iters++;

    #if DEBUG
        if (iters % 100 == 0) {
            cout << "iter " << iters << endl;
        }
    #endif

        // Select random (g, f) pair
        int g = uniform_int_distribution<int>(0, circt.G - 1)(rng);
        int f_idx = uniform_int_distribution<int>(0, circt.gates[g].p_faults.size() - 1)(rng);
        auto it = circt.gates[g].p_faults.begin();
        advance(it, f_idx);
        Fault f = it->first;

        Event E_prime = E;
        if (E[g] == Fault::NONE) {
            // Case 1: g not in any gate fault pair in E
            // -> E' = E U {(g, f)}
            E_prime[g] = f;
            double p_g = (circt.gates[g].p_failure > 0.0) ? circt.gates[g].p_failure : p;
            double q = min(1.0, (p_g / (1.0 - p_g)) * circt.gates[g].p_faults.at(f));
            if (RNG::bernoulli(q) && decode(circt, E_prime)) {
                E = E_prime;
            }
            // If trial failed or E' not in F, don't change E
        } else {
            // Case 2: there exists some h for which (g, h) in E
            // -> E' = (E U {(g, f)})\(g, h)
            int h = E[g];
            if (f == h) {
                // Acceptance probability (of removing fault) = 1, test if E' in F
                E_prime[g] = Fault::NONE;
                if (decode(circt, E_prime)) {
                    E = E_prime;
                }
                // If E' not in F, don't change E
            } else {
                // Replace (g, h) with (g, f) with probability P_g(f)
                double q = circt.gates[g].p_faults.at(f);
                if (RNG::bernoulli(q) && decode(circt, E_prime)) {
                    E = E_prime;
                }
                // If trial failed or E' not in F, don't change E
            }
        }

        // Check if current event E in F
        // if (decode(circt, E)) {
        //     failing_events.push_back(E);
        // }
        // We know E is in F
        failing_events.push_back(E);
    }

    auto end = chrono::high_resolution_clock::now();
    stats.metropolis_splits.push_back(chrono::duration_cast<chrono::milliseconds>(end - start));

    return failing_events;
}

// - Ratio Estimation -
static inline double g(double x) {
    return 1.0 / (1.0 + x);
}

double estimate_ratio(const Circuit& circt, const vector<Event>& failing_events, double p_i, double p_ip1) {
    vector<double> log_ratios;  // log(pi_i / pi_ip1) = -log(pi_ip1 / pi_p)
    for (const Event &E : failing_events) {
        log_ratios.push_back(log_event_probability_ratio(circt, E, p_i, p_ip1));
    }

    // Expectation computation helpers
    auto expectation_i = [&](double log_C) -> double {
        // mean(g(C * pi_i / pi_ip1))
        double sum = 0.0;
        for (double l_r : log_ratios) {
            sum += g(exp(log_C + l_r));
        }
        return sum / (double) log_ratios.size();
    };

    auto expectation_ip1 = [&](double log_C) -> double {
        // mean(g(C^-1 * pi_ip1 / pi_i))
        double sum = 0.0;
        for (double l_r : log_ratios) {
            sum += g(exp(-log_C - l_r));
        }
        return sum / (double) log_ratios.size();
    };

    // Solving for f(C) = Expectation_i(C) - Expectation_i+1(C) = 0
    auto f = [&](double log_C) -> double {
        return expectation_i(log_C) - expectation_ip1(log_C);
    };

    // Solve for C via root-finding
    double L = -60.0, R = 60.0;
    double fL = f(L), fR = f(R);
    // If signs are same, expand search range mildly or return best-guess ratio = exp((L+R)/2)
    if (fL * fR > 0) {
        // Try smaller bracket around 0
        L = -10; R = 10;
        fL = f(L); fR = f(R);
        if (fL * fR > 0) {
            // fallback: choose logC such that left ~= right approx via small grid search
            double best_logC = 0;
            double best_abs = fabs(f(0.0));
            for (double c=-20; c<=20; c+=0.5) {
                double val = fabs(f(c));
                if (val < best_abs) { best_abs = val; best_logC = c; }
            }
            cout << "Percentage error in approximation: " << (f(best_logC) / expectation_i(best_logC)) * 100  << "%" << endl;
            return std::exp(best_logC) * expectation_i(best_logC) / expectation_ip1(best_logC);
        }
    }

    // Bisection
    for (int iter=0; iter<200; ++iter) {
        double M = 0.5 * (L + R);
        double fM = f(M);
        if (fabs(fM) < 1e-12) {
            cout << "Percentage error in approximation: " << (f(M) / expectation_i(M)) * 100 << "%" << endl;
            return std::exp(M) * expectation_i(M) / expectation_ip1(M);
        }
        if (fL * fM <= 0) {
            R = M; fR = fM;
        } else {
            L = M; fL = fM;
        }
    }
    double logC = 0.5 * (L + R);
    cout << "Percentage error in approximation: " << (f(logC) / expectation_i(logC)) * 100 << "%" << endl;
    return std::exp(logC) * expectation_i(logC) / expectation_ip1(logC);
}

// --- Monte Carlo Simulation ---
double monte_carlo_initial_sim(int d, double p_phy) {
    auto start = chrono::high_resolution_clock::now();

    int shots = 10000000;

    ostringstream oss;
    oss << "./.venv/bin/python3 src/python/seed_monte_carlo.py " << d << " " << p_phy << " --shots " << shots;
    string cmd = oss.str();

    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        exit(1);
    }

    char buffer[128];
    string result = "";
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }

    pclose(pipe);

    double p_logical_0 = stod(result);
    cout << "Estimated logical error from Monte Carlo simulation: " << p_logical_0 << endl;

    auto end = chrono::high_resolution_clock::now();
    stats.mc_seed = chrono::duration_cast<chrono::milliseconds>(end - start);

    return p_logical_0;
}

// --- Rare Event Simulation ---
double rare_event_simulation(Circuit& circt, double p_0, double p_logical_target, int N) {
    auto start = chrono::high_resolution_clock::now();

    double p_logical_0 = monte_carlo_initial_sim(circt.d, p_0);

    // Physical error splitting sequence and estimated logical error rate vectors
    vector<double> p_seq;
    vector<double> p_logical;

    p_seq.push_back(p_0);
    p_logical.push_back(p_logical_0);

    // Get initial failing event to start MCMC
    // Event E_init = Event(circt.G, Fault::NONE);
    Event E_init = find_logical_error(circt);

    // Estimate ratio of each consecutive pair in sequence (p_i, p_i+1)
    int i = 0;
    while (p_logical.back() >= p_logical_target) {
    // while (p_seq.back() >= 1e-5) {
        double p_i = p_seq.back();

        // Get p_i+1 in splitting sequence
        double p_next = splitting_sequence_next(circt, p_i);
        p_seq.push_back(p_next);

        cout << "Level " << i << ": p_i = " << p_i << " -> p_i+1 = " << p_next << endl;

        // Run Metropolis routine to get N failing samples from pi_i|F
        vector<Event> failing_samples = run_metropolis_routine(circt, p_i, E_init, N);

        // Estimate ratio
        double ratio = estimate_ratio(circt, failing_samples, p_i, p_next);
        cout << "Estimated ratio: " << ratio << endl;

        // Estimate logical next failure in sequence with ratio
        double p_logical_est = p_logical.back() * ratio;
        p_logical.push_back(p_logical_est);
        cout << "p_logical_est[" << (i + 1) << "] = " << p_logical_est << endl;

        // Initial failing event for next level
        if (!failing_samples.empty()) {
            E_init = failing_samples.back();
        }

        i++;
    }

    cout << "Final estimate for target p = " << p_seq.back() << " is p_logical_est = " << p_logical.back() << endl;

    // Summarize (physical, logical) error rates
    cout << "--- Error Rates Summary ---" << endl;
    for (int i = 0; i < (int) p_seq.size(); ++i) {
        cout << "Level " << i << ": p_phy = " << p_seq[i] << ", p_logical = " << p_logical[i] << endl;
    }
    cout << "--- End Error Summary ---\n" << endl;

    auto end = chrono::high_resolution_clock::now();
    stats.total = chrono::duration_cast<chrono::milliseconds>(end - start);

    return p_logical.back();
}


// --- main() ---
// Arguments:
//   d: code distance
//   p_0: initial physical error rate
//   p_target: target logical error rate
//   N: per-level MCMC failures
int main(int argc, char* argv[]) {
    // Simulation parameters (defaults)
    int d = 3;
    double p = 0.0;
    double p_0 = 1e-3;
    double p_target = 1e-7;
    // int mc_required_failures = 1000;
    int N = 1000;

    // Parse arguments
    if (argc > 1) {
        d = stoi(argv[1]);
    }
    if (argc > 2) {
        p_0 = stod(argv[2]);
    }
    if (argc > 3) {
        p_target = stod(argv[3]);
    }
    if (argc > 4) {
        N = stoi(argv[4]);
    }

    Circuit circt = build_surface_code_circuit(d, p);

    double p_logical = rare_event_simulation(circt, p_0, p_target, N);

    cout << "Estimation of logical error rate: " << p_logical << endl;

    summarize_stats(stats);

    return 0;
}
