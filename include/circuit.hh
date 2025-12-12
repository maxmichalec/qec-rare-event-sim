#ifndef CIRCUIT_H_
#define CIRCUIT_H_

#pragma once

#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;

// --- Gate & Fault Model ---
enum GateType {
    INIT_X,
    INIT_Z,
    MEASURE_X,
    MEASURE_Z,
    CNOT,
    IDLE
};

enum Fault {
    NONE,
    X,
    Y,
    Z,
    XI,
    IX,
    ZI,
    IZ,
    XX,
    YY,
    ZZ
};

const map<Fault, string> faultToString = {
    {Fault::NONE, "NONE"},
    {Fault::X, "X"},
    {Fault::Y, "Y"},
    {Fault::Z, "Z"},
    {Fault::XI, "XI"},
    {Fault::IX, "IX"},
    {Fault::ZI, "ZI"},
    {Fault::IZ, "IZ"},
    {Fault::XX, "XX"},
    {Fault::YY, "YY"},
    {Fault::ZZ, "ZZ"}
};

struct Point {
    int x;
    int y;
};

Point operator+(Point p1, const Point& p2);

inline int point_to_idx(const Point& p, int d) {
    float y = p.y - fmod(p.x, 2);
    return (int) (p.x + y * (d + 0.5));
}

struct Gate {
    int id;
    int round = 0;
    GateType g_type = GateType::IDLE;
    Point qubit;
    int target = 0;

    double p_failure = 0.0;  // gate-unique failure probability
    // p_faults[k] = P(f = k | gate fail), sum(p_faults) = 1
    map<Fault, double> p_faults;
};

struct Circuit {
    int d = 3;
    int n_data_qubits = 0;
    int G = 0;  // gate count
    int gates_per_round = 0;
    vector<Gate> gates;
};

Circuit build_surface_code_circuit(int d, double p);

#endif /* CIRCUIT_H_ */
