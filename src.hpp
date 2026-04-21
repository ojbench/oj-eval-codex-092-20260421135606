#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"

#include <vector>

class resistive_network {
private:
    int interface_size;
    int connection_size;

    struct Edge { int u, v; fraction g; };
    std::vector<Edge> edges;
    // Laplacian matrix L (n x n)
    std::vector<std::vector<fraction>> L;

    // Solve (L without last row/col) * x = b for x (size n-1)
    std::vector<fraction> solve_reduced(const std::vector<fraction>& b) {
        int n = interface_size - 1;
        // Build A as top-left (n x n) submatrix of L
        std::vector<std::vector<fraction>> A(n, std::vector<fraction>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A[i][j] = L[i][j];
            }
        }
        std::vector<fraction> rhs = b;

        // Gaussian elimination
        for (int col = 0, row = 0; col < n && row < n; ++col, ++row) {
            int pivot = row;
            // Find non-zero pivot in or below current row
            while (pivot < n) {
                fraction zero(0);
                if (!(A[pivot][col] == zero)) break;
                ++pivot;
            }
            if (pivot == n) {
                // Should not happen for connected network with grounding
                throw resistive_network_error();
            }
            if (pivot != row) {
                std::swap(A[pivot], A[row]);
                std::swap(rhs[pivot], rhs[row]);
            }

            // Eliminate rows below
            fraction piv = A[row][col];
            for (int i = row + 1; i < n; ++i) {
                fraction zero(0);
                if (A[i][col] == zero) continue;
                fraction factor = A[i][col] / piv;
                for (int j = col; j < n; ++j) {
                    A[i][j] = A[i][j] - factor * A[row][j];
                }
                rhs[i] = rhs[i] - factor * rhs[row];
            }
        }

        // Back substitution
        std::vector<fraction> x(n, fraction(0));
        for (int i = n - 1; i >= 0; --i) {
            fraction sum(0);
            for (int j = i + 1; j < n; ++j) {
                sum = sum + A[i][j] * x[j];
            }
            fraction denom = A[i][i];
            fraction num = rhs[i] - sum;
            x[i] = num / denom;
        }
        return x;
    }

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        interface_size = interface_size_;
        connection_size = connection_size_;
        L.assign(interface_size, std::vector<fraction>(interface_size, fraction(0)));

        edges.reserve(connection_size);
        for (int k = 0; k < connection_size; ++k) {
            int u = from[k];
            int v = to[k];
            // Convert to 0-based indices
            int ui = u - 1;
            int vi = v - 1;
            fraction g = fraction(1) / resistance[k];
            edges.push_back({u, v, g});

            L[ui][ui] = L[ui][ui] + g;
            L[vi][vi] = L[vi][vi] + g;
            L[ui][vi] = L[ui][vi] - g;
            L[vi][ui] = L[vi][ui] - g;
        }
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int a = interface_id1;
        int b = interface_id2;
        int n = interface_size;

        // Build current vector for non-ground nodes (ground is node n)
        std::vector<fraction> rhs(n - 1, fraction(0));
        if (a != n) rhs[a - 1] = rhs[a - 1] + fraction(1);
        if (b != n) rhs[b - 1] = rhs[b - 1] - fraction(1);

        std::vector<fraction> volt = solve_reduced(rhs);

        fraction va = (a == n) ? fraction(0) : volt[a - 1];
        fraction vb = (b == n) ? fraction(0) : volt[b - 1];
        return va - vb; // With 1A injection, V difference equals equivalent resistance
    }

    fraction get_voltage(int id, fraction current[]) {
        int n = interface_size;
        // Build current vector for non-ground nodes
        std::vector<fraction> rhs(n - 1, fraction(0));
        for (int i = 0; i < n - 1; ++i) rhs[i] = current[i];
        std::vector<fraction> volt = solve_reduced(rhs);
        return volt[id - 1];
    }

    fraction get_power(fraction voltage[]) {
        fraction total(0);
        for (const auto &e : edges) {
            int ui = e.u - 1;
            int vi = e.v - 1;
            fraction dv = voltage[ui] - voltage[vi];
            total = total + e.g * (dv * dv);
        }
        return total;
    }
};


#endif //SRC_HPP

