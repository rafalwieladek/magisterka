/*
 * Kumari UAV Path Planner - Header File
 *
 * Pełna implementacja artykułu:
 * "Multi-UAV Path Planning for Connectivity-Based Sweep Coverage"
 * Kumari & Srirangarajan (2025), Ad Hoc Networks, 178, 103966
 *
 * Zmienne decyzyjne:
 * - z^n_{i,k}: Binarna, czy UAV k jest w punkcie siatki i w kroku n
 * - c^n_{i,k}: Ciągła, czy punkt siatki i jest pokryty przez UAV k w kroku n
 * - c_i: Ciągła, czy punkt siatki i jest pokryty przez jakikolwiek UAV
 *
 * Ograniczenia (równania 2-10):
 * (2) Pozycja UAV: Σ z^n_{i,k} = 1
 * (3) Kolizja: Σ z^n_{i,k} ≤ 1
 * (4) Pokrycie: Σ Σ z^n_{i,k} ≥ 1
 * (5) Mobilność: z^n_{i,k} ≤ Σ z^{n-1}_{j,k}
 * (6) Warunek początkowy
 * (7) Pokrycie punktu: c^n_{i,n} = Σ z^n_{j,n}
 * (8) Globalne pokrycie: c_i ≥ c^n_{i,n}
 * (9) Całkowite pokrycie: c_i ≤ Σ Σ c^n_{i,k}
 * (10) Ograniczenia zmiennych: 0 ≤ c_i, c^n_{i,k} ≤ 1
 */

#ifndef KUMARI_UAV_PATH_PLANNER_H
#define KUMARI_UAV_PATH_PLANNER_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;
typedef vector<vector<int>> IntMatrix;

// Makra do indeksowania zmiennych
// z^n_{i,k}: z[n][i][k] -> z_flat[n * num_grid * num_uavs + i * num_uavs + k]
#define Z_VAR(n, i, k, num_grid, num_uavs) ((n) * (num_grid) * (num_uavs) + (i) * (num_uavs) + (k))

// c^n_{i,k}: c_time[n][i][k]
#define C_TIME_VAR(n, i, k, num_grid, num_uavs) ((n) * (num_grid) * (num_uavs) + (i) * (num_uavs) + (k))

// c_i: c_global[i]
#define C_GLOBAL_VAR(i) (i)

class KumariUAVPathPlanner {
public:
    // Parametry problemu
    int grid_width;           // Szerokość siatki X
    int grid_height;          // Wysokość siatki Y
    int num_grid_points;      // Całkowita liczba punktów siatki (X * Y)
    int num_uavs;             // Liczba UAV (N)
    int max_time_steps;       // Maksymalna liczba kroków czasowych (K_max)

    double comm_range;        // Zasięg komunikacji (R_c)
    double sensing_range;     // Zasięg czujnika (R_s)
    double desired_coverage;  // Pożądany poziom pokrycia (cr)

    // Pozycje punktów siatki
    vector<pair<double, double>> grid_positions;  // (x, y) dla każdego punktu

    // Pozycja bazy
    pair<double, double> base_station;

    // Obszar zainteresowania
    vector<int> area_of_interest;  // Indeksy punktów w obszarze zainteresowania

    // CPLEX
    IloEnv env;
    IloModel model;
    IloCplex cplex;

    // Zmienne decyzyjne
    IloBoolVarArray z;           // z^n_{i,k} - pozycja UAV
    IloNumVarArray c_time;       // c^n_{i,k} - pokrycie w kroku
    IloNumVarArray c_global;     // c_i - globalne pokrycie

    // Pomiary
    struct Measurements {
        double build_time;
        double solve_time;
        double total_time;
        double objective_value;
        double mip_gap;
        int status;
        int nodes_processed;
        vector<vector<int>> uav_positions;  // Pozycje UAV w każdym kroku
        vector<int> coverage_per_step;      // Pokrycie w każdym kroku
    } measurements;

    // Konstruktor
    KumariUAVPathPlanner(int grid_width, int grid_height, int num_uavs,
        int max_time_steps, double comm_range, double sensing_range,
        double desired_coverage = 0.95);

    // Destruktor
    ~KumariUAVPathPlanner();

    // Metody publiczne
    void generate_grid_instance(int seed = 42);
    void set_base_station(double x, double y);
    void set_area_of_interest(const vector<int>& points);
    bool build_model();
    bool solve(double time_limit = 300.0);
    void print_results();
    void print_solution();

private:
    // Metody prywatne
    double calculate_distance(double x1, double y1, double x2, double y2);
    int get_grid_point(int x, int y);
    pair<int, int> get_grid_coordinates(int point_idx);

    // Budowanie ograniczeń
    void add_constraint_position_uav();           // Eq. (2)
    void add_constraint_collision();              // Eq. (3)
    void add_constraint_coverage();               // Eq. (4)
    void add_constraint_mobility();               // Eq. (5)
    void add_constraint_initial_condition();      // Eq. (6)
    void add_constraint_coverage_point();         // Eq. (7)
    void add_constraint_global_coverage();        // Eq. (8)
    void add_constraint_total_coverage();         // Eq. (9)
    void add_constraint_variable_bounds();        // Eq. (10)
};

#endif // KUMARI_UAV_PATH_PLANNER_H