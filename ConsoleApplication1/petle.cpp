#include "stdafx.h"
/*
 * ======================================================================================
 * SYSTEM: Multi-UAV SAR Mission Planner
 * MODEL: Pełne pętle z powrotem do bazy + Łączność
 * ======================================================================================
 */

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>

ILOSTLBEGIN

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

class SARMissionPlanner {
public:
    int N, K;
    double Tmax, R, V;

    Matrix pos;
    Vector priority, service_time, early_time, due_time;
    Matrix dist_matrix;

    IloEnv env;
    IloModel model;
    IloCplex cplex;

    // x[k][i][j] - pełna trasa (z bazy, przez cele, z powrotem do bazy)
    IloArray<IloBoolVarArray> x;
    // y[k][i] - czy dron k "zalicza" węzeł i
    IloArray<IloBoolVarArray> y;
    // t[k][i] - czas przybycia do węzła 
    IloArray<IloNumVarArray> t;

    // Zmienne łączności (muszą być spełnione, gdy drony "mijają" punkty)
    IloArray<IloBoolVarArray> r;
    IloArray<IloNumVarArray> f;

    SARMissionPlanner(int n, int k, double tmax, double radius, double velocity)
        : N(n), K(k), Tmax(tmax), R(radius), V(velocity),
        env(), model(env), cplex(model)
    {
        pos.assign(N, Vector(2));
        priority.assign(N, 0.0);
        service_time.assign(N, 0.0);
        early_time.assign(N, 0.0);
        due_time.assign(N, Tmax);
        dist_matrix.assign(N, Vector(N));

        x = IloArray<IloBoolVarArray>(env, K);
        y = IloArray<IloBoolVarArray>(env, K);
        t = IloArray<IloNumVarArray>(env, K);
        r = IloArray<IloBoolVarArray>(env, K);
        f = IloArray<IloNumVarArray>(env, N);

        for (int k = 0; k < K; k++) {
            x[k] = IloBoolVarArray(env, N * N);
            y[k] = IloBoolVarArray(env, N);
            t[k] = IloNumVarArray(env, N, 0.0, Tmax);
            r[k] = IloBoolVarArray(env, N * N);
        }
        for (int i = 0; i < N; i++) {
            f[i] = IloNumVarArray(env, N, 0.0, (double)(N * K));
        }
    }

    ~SARMissionPlanner() { env.end(); }

    double calc_dist(int i, int j) const {
        return std::sqrt(std::pow(pos[i][0] - pos[j][0], 2) + std::pow(pos[i][1] - pos[j][1], 2));
    }

    void generate_complex_scenario() {
        srand(999);
        pos[0] = { 10, 10 };
        priority[0] = 0; service_time[0] = 0;
        early_time[0] = 0.0; due_time[0] = Tmax;

        // Cel Strategiczny 
        pos[1] = { 90, 90 };
        priority[1] = 500;
        service_time[1] = 10.0;
        early_time[1] = 30.0;  // Wąskie okno
        due_time[1] = 35.0;

        // Drugi cel
        pos[2] = { 50, 60 };
        priority[2] = 200;
        service_time[2] = 5.0;
        early_time[2] = 15.0;  // Wąskie okno
        due_time[2] = 25.0;

        // Cele poboczne / Potencjalne węzły przekaźnikowe
        for (int i = 3; i < N; i++) {
            pos[i][0] = rand() % 80 + 10;
            pos[i][1] = rand() % 80 + 10;
            priority[i] = 10 + rand() % 20;
            service_time[i] = 2.0;
            early_time[i] = 0.0;
            due_time[i] = Tmax; // Luźne okna
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                dist_matrix[i][j] = calc_dist(i, j);
            }
        }
    }

    void build_model() {
        // [1] FUNKCJA CELU
        IloExpr obj(env);
        for (int k = 0; k < K; k++) {
            for (int i = 1; i < N; i++) {
                obj += priority[i] * y[k][i];
            }
        }
        model.add(IloMaximize(env, obj));

        // [2] ROUTING (Pełne Pętle )
        for (int k = 0; k < K; k++) {
            // Wylot z bazy
            IloExpr start(env);
            for (int j = 1; j < N; j++) start += x[k][0 * N + j];
            model.add(start <= 1); // Dron nie musi startować

            // Powrót do bazy
            IloExpr stops(env);
            for (int i = 1; i < N; i++) stops += x[k][i * N + 0];

            // Bilans bazy: Zawsze musi wrócić, jeśli wystartował (Closed Loop!)
            model.add(stops == start);

            // Bilans dla węzłów pośrednich: WLOT = WYLOT
            // To jest główna różnica! Dron wpada i wylatuje (nie zatrzymuje się ostatecznie).
            for (int i = 1; i < N; i++) {
                IloExpr flow_in(env);
                IloExpr flow_out(env);
                for (int j = 0; j < N; j++) {
                    if (i != j) {
                        flow_in += x[k][j * N + i];
                        flow_out += x[k][i * N + j];
                    }
                }
                // Dron musi wylecieć, jeśli wleciał
                model.add(flow_in == flow_out);

                // Powiązanie zmiennej 'y' (czy zaliczył cel)
                // Jeśli przez niego przeleciał, to go zaliczył.
                // Pozwalamy na wielokrotne przeloty (flow_out >= y) - to pozwala wrócić po własnych śladach.
                model.add(flow_out >= y[k][i]);
                model.add(flow_out <= (double)N * y[k][i]);
            }

            // Baza nie jest liczona jako cel
            model.add(y[k][0] == 0);
        }

        // Limit kolizji: Punkt może być zaliczony przez max 1 drona
        // Zakomentowano, aby umożliwić lot kilku dronów przez ten sam punkt "gęsiego".
         /*for(int i=1; i<N; i++) {
             IloExpr occ(env); 
             for(int k=0; k<K; k++) occ += y[k][i];
             model.add(occ <= 1);
         }*/

        // [3] CZAS I OKNA CZASOWE 
        double M_time = Tmax + 100.0;

        for (int k = 0; k < K; k++) {
            model.add(t[k][0] == 0); // Start

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (i != j) {
                        //zmiana
                        if (j != 0) {
                            double travel = dist_matrix[i][j] / V;
                            double s_time = (i == 0) ? 0 : service_time[i];
                            // Czas rośnie wraz z lotem
                            model.add(t[k][j] >= t[k][i] + travel + s_time - M_time * (1 - x[k][i * N + j]));
                        }
                    }
                }

                // Okna czasowe aktywne, jeśli cel został zaliczony (y=1)
                if (i > 0) {
                    model.add(t[k][i] >= early_time[i] * y[k][i]);
                    model.add(t[k][i] <= due_time[i] + M_time * (1 - y[k][i]));
                }

                // Limit misji
                model.add(t[k][i] <= Tmax);
            }
        }

        // [4] ŁĄCZNOŚĆ (Kumari)
        // Ograniczenie dystansu
        for (int k = 0; k < K; k++) for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
            if (dist_matrix[i][j] > R)
                model.add(r[k][i * N + j] == 0);

        // Łącze radiowe r[i][j] istnieje, jeśli przez i oraz j "przeleciały" drony.
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                IloExpr rs(env);
                for (int k = 0; k < K; k++) rs += r[k][i * N + j];

                IloExpr yi(env), yj(env);
                if (i == 0) yi += 1; else for (int k = 0; k < K; k++) yi += y[k][i];
                if (j == 0) yj += 1; else for (int k = 0; k < K; k++) yj += y[k][j];

                model.add(rs <= yi);
                model.add(rs <= yj);
            }
        }

        // [5] WIRTUALNY PRZEPŁYW DANYCH (Commodity Flow)
        IloExpr source_out(env);
        for (int j = 1; j < N; j++) {
            source_out += f[0][j];
            source_out -= f[j][0];
        }

        IloExpr total_demand(env);
        for (int i = 1; i < N; i++) for (int k = 0; k < K; k++) total_demand += y[k][i];

        model.add(source_out == total_demand);

        for (int i = 1; i < N; i++) {
            IloExpr flow_net(env);
            for (int j = 0; j < N; j++) {
                flow_net += f[j][i];
                flow_net -= f[i][j];
            }
            IloExpr visited(env);
            for (int k = 0; k < K; k++) visited += y[k][i];

            model.add(flow_net == visited);
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                IloExpr r_total(env);
                for (int k = 0; k < K; k++) r_total += r[k][i * N + j];
                model.add(f[i][j] <= (double)(N * K) * r_total);
            }
        }

        // Zmuszenie do zaliczenia Celu strategicznego (Węzeł 1)
        IloExpr super_target(env);
        for (int k = 0; k < K; k++) super_target += y[k][1];
        model.add(super_target >= 1);
    }

    void solve() {
        //cplex.setOut(env.getNullStream());
        cplex.setParam(IloCplex::TiLim, 240);
        cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);

        std::cout << "\n[INFO] Rozpoczynam optymalizacje ..." << std::endl;

        if (cplex.solve()) {
            std::cout << ">>> SUKCES! ZNALEZIONO ROZWIAZANIE." << std::endl;
            std::cout << "Wartosc funkcji celu: " << cplex.getObjValue() << std::endl;
            printRoutes();
            visualize();
        }
        else {
            std::cout << "[BLAD] INFEASIBLE." << std::endl;
        }
    }

    void printRoutes() {
        std::cout << "\n=== SZCZEGOLY TRAS (Weryfikacja powrotu) ===\n";
        for (int k = 0; k < K; k++) {
            bool active = false;
            for (int j = 1; j < N; j++) if (cplex.getValue(x[k][0 * N + j]) > 0.5) active = true;

            if (!active) {
                std::cout << "Dron " << k << ": [BAZA - nieaktywny]\n";
                continue;
            }

            std::cout << "Dron " << k << ": 0";
            int curr = 0;
            int steps = 0;

            while (steps < N * 3) { // Zabezpieczenie pętli
                int next = -1;
                for (int j = 0; j < N; j++) {
                    if (curr == j) continue;
                    try { if (cplex.getValue(x[k][curr * N + j]) > 0.5) { next = j; break; } }
                    catch (...) {}
                }

                if (next != -1) {
                    std::cout << " -> " << next;
                    if (next == 0) {
                        std::cout << " [LADOWANIE]";
                        break; // Wrócił do bazy = koniec trasy
                    }
                    curr = next;
                    steps++;
                }
                else {
                    break;
                }
            }
            std::cout << "\n";
        }
        std::cout << "===========================================\n";
    }

    void visualize() {
        std::ofstream file("sar_closed_loop.svg");
        if (!file) return;

        double s = 6.0;
        int w = 800, h = 800;
        int off = 50;

        file << "<svg width=\"" << w << "\" height=\"" << h << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"#fdfdfd\" />\n";

        // Łącza radiowe (Backbone)
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double val = 0;
                try { for (int k = 0; k < K; k++) val += cplex.getValue(r[k][i * N + j]); }
                catch (...) {}
                if (val > 0.5) {
                    file << "<line x1=\"" << (off + pos[i][0] * s) << "\" y1=\"" << (h - off - pos[i][1] * s)
                        << "\" x2=\"" << (off + pos[j][0] * s) << "\" y2=\"" << (h - off - pos[j][1] * s)
                        << "\" stroke=\"#b3b3b3\" stroke-width=\"3\" stroke-dasharray=\"6,6\"/>\n";
                }
            }
        }

        // Trasy dronów
        std::string cols[] = { "red", "green", "blue", "orange", "purple", "brown", "yellow", "pink"};
        for (int k = 0; k < K; k++) {
            double dx = (k % 2 == 0 ? 1 : -1) * (k / 2 + 1) * 3.0;
            double dy = (k % 2 == 0 ? -1 : 1) * (k / 2 + 1) * 3.0;

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    double v = 0;
                    try { v = cplex.getValue(x[k][i * N + j]); }
                    catch (...) {}

                    if (i != j && v > 0.5) {
                        double x1 = off + pos[i][0] * s + dx;
                        double y1 = h - off - pos[i][1] * s + dy;
                        double x2 = off + pos[j][0] * s + dx;
                        double y2 = h - off - pos[j][1] * s + dy;

                        // Strzałka na 80% długości by pokazać kierunek powrotny
                        double mx = x1 + 0.8 * (x2 - x1);
                        double my = y1 + 0.8 * (y2 - y1);

                        file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2
                            << "\" stroke=\"" << cols[k % 6] << "\" stroke-width=\"2\"/>\n";

                        file << "<circle cx=\"" << mx << "\" cy=\"" << my << "\" r=\"3\" fill=\"" << cols[k % 6] << "\"/>\n";
                    }
                }
            }
        }

        for (int i = 0; i < N; i++) {
            double cx = off + pos[i][0] * s;
            double cy = h - off - pos[i][1] * s;

            std::string fill = (i == 1) ? "gold" : (i == 0 ? "red" : "black");
            double r_size = (i == 1) ? 10 : 5;

            if (i == 0) {
                file << "<rect x=\"" << cx - 10 << "\" y=\"" << cy - 10 << "\" width=\"20\" height=\"20\" fill=\"red\"/>";
            }
            else {
                file << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r_size << "\" fill=\"" << fill << "\" stroke=\"black\"/>";
            }
            file << "<text x=\"" << cx + 10 << "\" y=\"" << cy << "\" font-size=\"11\" font-weight=\"bold\">Node " << i << "</text>\n";
        }

        file << "</svg>\n";
        file.close();

        std::cout << "[INFO] Zapisano plik: sar_closed_loop.svg" << std::endl;
        system("start sar_closed_loop.svg");
    }
};

int main() {
    try {
        SARMissionPlanner planner(20, 3, 60.0, 35.0, 10.0);
        planner.generate_complex_scenario();
        planner.build_model();
        planner.solve();
    }
    catch (IloException& e) { std::cerr << "CPLEX ERROR: " << e << "\n"; }
    return 0;
}