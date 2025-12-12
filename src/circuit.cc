#include "circuit.hh"

Point operator+(Point p1, const Point& p2) {
    Point p_sum;
    p_sum.x = p1.x + p2.x;
    p_sum.y = p1.y + p2.y;
    return p_sum;
}

Circuit build_surface_code_circuit(int d, double p) {
    Circuit circt;
    circt.d = d;
    circt.n_data_qubits = d * d;
    circt.G = 0;
    circt.gates_per_round = 0;

    vector<Gate> gates;
    // int n_ancilla = (d * d) - 1;

    vector<pair<Point, bool>> ancilla;
    for (int x = 0; x <= d; ++x) {
        for (int y = 0; y <= d; ++y) {
            Point q = { .x = x * 2, .y = y * 2 };
            bool on_boundary_1 = (x == 0) || (x == d);
            bool on_boundary_2 = (y == 0) || (y == d);
            bool parity = (x % 2) != (y % 2);
            if (on_boundary_1 && parity) {
                continue;
            }
            if (on_boundary_2 && !parity) {
                continue;
            }
            // if parity: MEASURE_X else: MEASURE_Z
            ancilla.push_back(pair<Point, bool>(q, parity));
        }
    }

    // Enumerate data qubit ids to ensure CNOTs are valid
    vector<int> data_qubits;
    for (float x = 0.5; x < d; x += 1.0) {
        for (float y = 0.5; y < d; y += 1.0) {
            Point q = { .x = (int) (x * 2.0), .y = (int) (y * 2.0) };
            int q_id = point_to_idx(q, d);
            data_qubits.push_back(q_id);
        }
    }

    vector<Point> x_order = {
        { .x = 1, .y = 1 },
        { .x = -1, .y = 1 },
        { .x = 1, .y = -1 },
        { .x = -1, .y = -1 }
    };
    vector<Point> z_order = {
        { .x = 1, .y = 1 },
        { .x = 1, .y = -1 },
        { .x = -1, .y = 1 },
        { .x = -1, .y = -1 }
    };

    // For each round
    for (int r = 0; r < d; ++r) {
        // Prep
        for (pair<Point, bool>& a : ancilla) {
            Gate gate;
            gate.round = r;
            gate.qubit = a.first;
            gate.p_failure = p;
            if (a.second) {
                gate.g_type = GateType::INIT_X;
            } else {
                gate.g_type = GateType::INIT_Z;
            }
            gate.p_faults.insert(pair<Fault, double>(Fault::X, 1.0 / 3.0));
            gate.p_faults.insert(pair<Fault, double>(Fault::Y, 1.0 / 3.0));
            gate.p_faults.insert(pair<Fault, double>(Fault::Z, 1.0 / 3.0));
            gates.push_back(gate);
        }

        // CNOTs
        for (pair<Point, bool>& a : ancilla) {
            // 4 CNOTs per stabalizer
            for (int i = 0; i < 4; ++i) {
                Point data;
                if (a.second) { // X
                    data = a.first + x_order[i];
                } else { // Z
                    data = a.first + z_order[i];
                }

                int data_id = point_to_idx(data, d);
                bool data_found = false;
                for (int q_id : data_qubits) {
                    if (data_id == q_id) {
                        data_found = true;
                        break;
                    }
                }
                if (data_found) {
                    Gate gate;
                    gate.round = r;
                    gate.qubit = a.first;
                    gate.target = i;
                    gate.p_failure = p;
                    gate.g_type = GateType::CNOT;
                    gate.p_faults.insert(pair<Fault, double>(Fault::XI, 1.0 / 7.0));
                    gate.p_faults.insert(pair<Fault, double>(Fault::IX, 1.0 / 7.0));
                    gate.p_faults.insert(pair<Fault, double>(Fault::ZI, 1.0 / 7.0));
                    gate.p_faults.insert(pair<Fault, double>(Fault::IZ, 1.0 / 7.0));
                    gate.p_faults.insert(pair<Fault, double>(Fault::XX, 1.0 / 7.0));
                    gate.p_faults.insert(pair<Fault, double>(Fault::YY, 1.0 / 7.0));
                    gate.p_faults.insert(pair<Fault, double>(Fault::ZZ, 1.0 / 7.0));
                    gates.push_back(gate);
                }
                // else: Edge ancilla, does not have CNOT for round i
            }
        }

        // Measure
        for (pair<Point, bool>& a : ancilla) {
            Gate gate;
            gate.round = r;
            gate.qubit = a.first;
            gate.p_failure = p;
            if (a.second) {
                gate.g_type = GateType::MEASURE_X;
                gate.p_faults.insert(pair<Fault, double>(Fault::Z, 1.0 / 2.0));
            } else {
                gate.g_type = GateType::MEASURE_Z;
                gate.p_faults.insert(pair<Fault, double>(Fault::X, 1.0 / 2.0));
            }
            gate.p_faults.insert(pair<Fault, double>(Fault::Y, 1.0 / 2.0));
            gates.push_back(gate);
        }
    }

    circt.gates = gates;
    circt.G = gates.size();
    circt.gates_per_round = gates.size() / d;

    return circt;
}
