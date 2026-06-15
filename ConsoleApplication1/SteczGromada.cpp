#include "stdafx.h"
/*
 * ======================================================================================
 * SYSTEM: Multi-UAV SAR Mission Planner
 * MODEL: Hybrydowy (Stecz & Gromada + Kumari & Srirangarajan)
 * ======================================================================================
 */

#include <ilcplex/ilocplex.h>   // Główna biblioteka IBM ILOG CPLEX (Concert Technology)
#include <iostream>             // Obsługa strumieni wejścia/wyjścia (konsola)
#include <vector>               // Dynamiczne tablice STD
#include <cmath>                // Funkcje matematyczne (np. std::sqrt, std::pow do odległości)
#include <fstream>              // Obsługa operacji na plikach (do zapisu pliku .svg)
#include <string>               // Obsługa ciągów znaków
#include <sstream>              // Konwersja typów (np. wczytywanie danych z klawiatury)
#include <iomanip>              // Formator dla liczb zmiennoprzecinkowych (std::setprecision)
#include <cstdlib>              // Standardowe funkcje systemowe (np. system("start file.svg"))

ILOSTLBEGIN // Makro CPLEX: definiuje niezbędne przestrzenie nazw (m.in. std i ilog)

// Aliasy dla tablic jednowymiarowych (Vector) i dwuwymiarowych (Matrix) typu double.

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

class SARMissionPlanner {
public:
    // =========================================================================
    // PARAMETRY SYSTEMOWE (DANE WEJŚCIOWE)
    // =========================================================================
    int N;          // Liczba węzłów na mapie (węzeł 0 = stacja bazowa, 1..N-1 = cele/przekaźniki)
    int K;          // Liczba dronów we flocie
    double Tmax;    // Budżet czasu/pojemność baterii (maksymalny czas trwania misji dla drona)
    double R;       // Promień (zasięg) skutecznej komunikacji radiowej między dronami
    double V;       // Prędkość przelotowa drona

    // =========================================================================
    // STRUKTURY DANYCH MAPY I CELÓW
    // =========================================================================
    Matrix pos;           // Tablica dwuwymiarowa [N][2] przechowująca współrzędne (X, Y) węzłów
    Vector priority;      // Wektor priorytetów p_i celów (im wyższa liczba, tym cel ważniejszy)
    Vector service_time;  // Czas potrzebny na zeskanowanie celu radarem SAR w węźle i
    Vector early_time;    // Początek okna czasowego (Earliest time) - przed tym czasem dron musi czekać
    Vector due_time;      // Koniec okna czasowego (Due date) - po tym czasie wykonanie zadania jest nieważne
    Matrix dist_matrix;   // Wstępnie obliczona macierz odległości między każdym węzłem i oraz j

    // =========================================================================
    // OBIEKTY SILNIKA CPLEX
    // =========================================================================
    IloEnv env;       // Środowisko CPLEX - zarządza pamięcią wszystkich zmiennych i ograniczeń
    IloModel model;   // Model matematyczny - kontener do którego wkładamy funkcję celu i ograniczenia
    IloCplex cplex;   // Sam solver - algorytm, który rozwiązuje model (znajduje minimum/maksimum)

    // =========================================================================
    // ZMIENNE DECYZYJNE MILP
    // =========================================================================
    // IloBoolVarArray to wektor zmiennych przyjmujących tylko 0 lub 1 (prawda/fałsz).
    // IloNumVarArray to wektor zmiennych ciągłych (np. liczby rzeczywiste).

    IloArray<IloBoolVarArray> x; // Trasa przelotu. x[k][i][j] = 1, jeśli dron k leci z i do j.
    IloArray<IloBoolVarArray> y; // Zawis docelowy. y[k][i] = 1, jeśli dron k zatrzymuje się w węźle i (sensor/przekaźnik).
    IloArray<IloNumVarArray> t;  // Harmonogram. t[k][i] - moment (sekunda) przybycia drona k do węzła i.

    IloArray<IloBoolVarArray> r; // Topologia sieci. r[k][i][j] = 1, jeśli zestawiono aktywne łącze radiowe między i, j.
    IloArray<IloNumVarArray> f;  // Przepływ wirtualny. f[i][j] - wielkość danych przesyłanych od i do j.

    // -------------------------------------------------------------------------
    // KONSTRUKTOR KLASY
    // Inicjalizuje rozmiary tablic i rejestruje zmienne w środowisku CPLEX.
    // -------------------------------------------------------------------------
    SARMissionPlanner(int n, int k, double tmax, double radius, double velocity)
        : N(n), K(k), Tmax(tmax), R(radius), V(velocity),
        env(), model(env), cplex(model)
    {
        // Alokacja pamięci dla struktur wejściowych
        pos.assign(N, Vector(2));
        priority.assign(N, 0.0);
        service_time.assign(N, 0.0);
        early_time.assign(N, 0.0);
        due_time.assign(N, Tmax);
        dist_matrix.assign(N, Vector(N));

        // Alokacja pamięci dla zmiennych w CPLEX.
        x = IloArray<IloBoolVarArray>(env, K);
        y = IloArray<IloBoolVarArray>(env, K);
        t = IloArray<IloNumVarArray>(env, K);
        r = IloArray<IloBoolVarArray>(env, K);
        f = IloArray<IloNumVarArray>(env, N);

        // Inicjalizacja zmiennych specyficznych dla drona (indeksowanych k)
        for (int k = 0; k < K; k++) {
            // Rejestracja w CPLEX zmiennych binarnych (N*N krawędzi, N wierzchołków)
            x[k] = IloBoolVarArray(env, N * N);
            y[k] = IloBoolVarArray(env, N);
            // Rejestracja zmiennych ciągłych (czasu). Granice to 0 oraz Tmax.
            t[k] = IloNumVarArray(env, N, 0.0, Tmax);
            r[k] = IloBoolVarArray(env, N * N);
        }

        // Zmienna przepływu jest wspólna dla całej sieci.
        // Górna granica przepływu to N*K (maksymalna teoretyczna liczba strumieni).
        for (int i = 0; i < N; i++) {
            f[i] = IloNumVarArray(env, N, 0.0, (double)(N * K));
        }
    }

    // -------------------------------------------------------------------------
    // DESTRUKTOR KLASY
    // -------------------------------------------------------------------------
    ~SARMissionPlanner() {
        env.end();
    }

    // -------------------------------------------------------------------------
    // FUNKCJA POMOCNICZA: Obliczanie odległości euklidesowej z twierdzenia Pitagorasa.
    // -------------------------------------------------------------------------
    double calc_dist(int i, int j) const {
        return std::sqrt(std::pow(pos[i][0] - pos[j][0], 2) + std::pow(pos[i][1] - pos[j][1], 2));
    }

    // -------------------------------------------------------------------------
    // GENERATOR DANYCH: Scenariusz "Głębokie Rozpoznanie" (Deep Recon)
    // -------------------------------------------------------------------------
    void generate_complex_scenario() {
        srand(999); // Ustawienie sztywnego ziarna dla uzyskania powtarzalnych wyników

        // Punkt (0) - Stacja Bazowa (SINK)
        pos[0] = { 10, 10 };
        priority[0] = 0; service_time[0] = 0;

        // Punkt (1) - Cel HVT (High Value Target). Oddalony i wysoko punktowany (p=500).
        // Wymusza na solverze zbudowanie długiego łańcucha przekaźników.
        pos[1] = { 90, 90 };
        priority[1] = 500;
        service_time[1] = 10.0;
        early_time[1] = 20.0; // Cel HVT jest niedostępny przez pierwsze 200 sekund (np. ukryty za przeszkodą)
        due_time[1] = 40.0;  // Po 40 sekundach cel HVT ucieka (np. przemieszcza się poza zasięg)

        // Punkty (2..N-1) - Cele poboczne i punkty przekaźnikowe. Niska punktacja.
        for (int i = 2; i < N; i++) {
            pos[i][0] = rand() % 80 + 10;
            pos[i][1] = rand() % 80 + 10;
            priority[i] = 10 + rand() % 20;
            service_time[i] = 2.0;
        }
        early_time[2] = 10.0;
        due_time[2] = 45.0;
        // Zbudowanie macierzy odległości dla szybszego dostępu w modelu
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                dist_matrix[i][j] = calc_dist(i, j);
            }
        }
    }

    // -------------------------------------------------------------------------
    // GŁÓWNY MODUŁ MATEMATYCZNY (Budowa ograniczeń MILP)
    // -------------------------------------------------------------------------
    void build_model() {

        // =========================================================
        // 1. FUNKCJA CELU 
        // Cel: Zmaksymalizować sumę priorytetów zdobytych celów.
        // =========================================================
        IloExpr obj(env);
        for (int k = 0; k < K; k++) {
            for (int i = 1; i < N; i++) {
                // Dodaje do sumy: (priorytet) razy (zmienna binarna 1/0)
                obj += priority[i] * y[k][i];
            }
        }
        model.add(IloMaximize(env, obj)); // Ustawia funkcję jako maksymalizację
        obj.end(); // Zwalnia pamięć tymczasowego wyrażenia

        // =========================================================
        // 2. WARSTWA ROUTINGU wg Stecz i Gromada
        // =========================================================
        for (int k = 0; k < K; k++) {
            // A. Ograniczenie Wylotu z Bazy (Węzeł 0)
            IloExpr start(env);
            for (int j = 1; j < N; j++) start += x[k][0 * N + j];

            // Dron NIE MUSI startować z bazy (relaksacja z '==' na '<=').
            model.add(start <= 1);

            // B. Gwarancja docelowej pozycji zawisu (Lądowania)
            IloExpr stops(env);
            for (int i = 1; i < N; i++) stops += y[k][i];

            // Jeżeli wyleciał z bazy (start=1), to musi zająć dokładnie jeden węzeł (stops=1)
            model.add(stops == start);

            // C. Ciągłość grafu
            for (int i = 1; i < N; i++) {
                IloExpr flow(env);
                for (int j = 0; j < N; j++) flow += x[k][j * N + i]; // Wlot drona do węzła i
                for (int j = 0; j < N; j++) flow -= x[k][i * N + j]; // Wylot drona z węzła i

                // Bilans = 1, jeśli dron stacjonuje w węźle i (y=1). 
                // Bilans = 0, jeśli tylko przelatuje (Wlot = Wylot).
                model.add(flow == y[k][i]);
                flow.end();
            }

            // Baza jako punkt docelowy po wylocie jest zabroniona
            model.add(y[k][0] == 0);
            start.end(); stops.end();
        }

        // Limit kolizji: dany punkt w przestrzeni powietrznej zajmuje maks 1 dron na raz
        for (int i = 1; i < N; i++) {
            IloExpr occ(env);
            for (int k = 0; k < K; k++) occ += y[k][i];
            model.add(occ <= 1);
            occ.end();
        }

        // =========================================================
        // 3. WARSTWA CZASOWA I OKNA
        // =========================================================
        // Dynamiczne wyliczenie dużej stałej M
        double M_time = Tmax + 100.0;

        for (int k = 0; k < K; k++) {
            model.add(t[k][0] == 0); // Każdy dron startuje w chwili t=0

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (i != j) {
                        double travel = dist_matrix[i][j] / V;
                        model.add(t[k][j] >= t[k][i] + travel - M_time * (1 - x[k][i * N + j]));
                    }
                }

                // Budżet baterii (rezerwacja energii na powrót).
                // Koszt to (Czas dolotu * 2) + (Czas pracy z radarem SAR).
                IloExpr cost(env);
                cost += t[k][i] * 2.0;
                cost += service_time[i];
                model.add(cost <= Tmax + M_time * (1 - y[k][i]));
                cost.end();

                // Okna Czasowe (Time Windows).
                if (i > 0) {
                    model.add(t[k][i] >= early_time[i] * y[k][i]); // Czekaj na otwarcie
                    model.add(t[k][i] <= due_time[i] + M_time * (1 - y[k][i])); // Nie spóźnij się
                }
            }
        }

        // =========================================================
        // 4. WARSTWA SIECIOWA (Łączność wg Kumari)
        // =========================================================
        // Fizyka propagacji fali. Jeżeli dystans > Promień R, łącze radiowe nie istnieje.
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (dist_matrix[i][j] > R) {
                        model.add(r[k][i * N + j] == 0);
                    }
                }
            }
        }

        // Łącze radiowe istnieje tylko wtedy,
        // gdy na OBU końcach (węzeł i, węzeł j) stacjonują drony.
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) continue;

                IloExpr rs(env); // Sumaryczna liczba aktywnych łączy i-j
                for (int k = 0; k < K; k++) rs += r[k][i * N + j];

                IloExpr yi(env), yj(env); // Sumaryczna liczba dronów w węźle i oraz j

                if (i == 0) yi += 1; else for (int k = 0; k < K; k++) yi += y[k][i]; // Baza aktywna zawsze
                if (j == 0) yj += 1; else for (int k = 0; k < K; k++) yj += y[k][j];

                // Liczba łączy (rs) nie może być większa niż liczba dronów w punkcie.
                // Jeżeli węzeł jest pusty (y=0), wtedy zamknięcie łącza (rs=0).
                model.add(rs <= yi);
                model.add(rs <= yj);

                rs.end(); yi.end(); yj.end();
            }
        }

        // =========================================================
        // 5. WARSTWA PRZEPŁYWU DANYCH (Commodity Flow)
        // =========================================================
        // Baza jako źródło sieci produkujące sygnał
        IloExpr source_out(env);
        for (int j = 1; j < N; j++) {
            source_out += f[0][j]; // Wypływ z bazy w stronę innych węzłów
            source_out -= f[j][0]; // Wpływ zwrotny (nie powinno go być)
        }

        // Baza musi wyprodukować tyle danych, ilu jest "odbiorców" (dronów poza bazą)
        IloExpr total_demand(env);
        for (int i = 1; i < N; i++) {
            for (int k = 0; k < K; k++) total_demand += y[k][i];
        }
        model.add(source_out == total_demand);
        source_out.end(); total_demand.end();

        // Bilans dla węzłów docelowych w terenie 
        for (int i = 1; i < N; i++) {
            IloExpr flow_net(env);
            for (int j = 0; j < N; j++) {
                flow_net += f[j][i]; // Sygnał wpływający do węzła
                flow_net -= f[i][j]; // Sygnał wypływający dalej
            }
            IloExpr visited(env);
            for (int k = 0; k < K; k++) visited += y[k][i];

            // Węzeł pochłania 1 jednostkę sygnału, jeżeli stoi w nim dron. 
            // Inaczej (visited=0), bilans wychodzi na zero (dane tylko przelatują).
            model.add(flow_net == visited);
            flow_net.end(); visited.end();
        }

        // Ograniczenie przepustowości
        // Sygnał (f) może płynąć z i do j tylko wtedy, gdy uruchomiono tam łącze radiowe (r).
        double M_flow = (double)(N * K);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                IloExpr r_total(env);
                for (int k = 0; k < K; k++) r_total += r[k][i * N + j];

                model.add(f[i][j] <= M_flow * r_total);
                r_total.end();
            }
        }
    }

    // -------------------------------------------------------------------------
    // ROZWIĄZYWANIE PROBLEMU (Silnik CPLEX)
    // -------------------------------------------------------------------------
    void solve() {
        //cplex.setOut(env.getNullStream()); // Wyciszenie potoku informacji technicznych z solvera
        cplex.setParam(IloCplex::TiLim, 120); // Limit czasu 120 sekund na instancję

        // Wyłączenie Pre-solve.
        // Zapobiega to usuwaniu przez optymalizator "nieużywanych" zmiennych z pamięci
        // i chroni przed błędami Memory Access
        cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);

        std::cout << "\n[INFO] Rozpoczynam optymalizacje CPLEX (Metoda Branch-and-Cut)..." << std::endl;

        if (cplex.solve()) {
            std::cout << ">>> SUKCES! ZNALEZIONO ROZWIAZANIE." << std::endl;
            std::cout << "Wartosc funkcji celu (Zdobyte punkty): " << cplex.getObjValue() << std::endl;
            visualize();
        }
        else {
            std::cout << "[BLAD] INFEASIBLE. Problem matematycznie sprzeczny." << std::endl;
        }
    }

    // -------------------------------------------------------------------------
    // MODUŁ WIZUALIZACJI: Generowanie pliku wektorowego SVG
    // -------------------------------------------------------------------------
    void visualize() {
        std::ofstream file("sar_mission_map.svg");
        if (!file) return;

        // Ustawienia przestrzeni rysunkowej
        double s = 6.0;         // Skala - powiększa obszar 100x100 m do 600x600 px
        int w = 800, h = 800;   // Rozmiar obrazu w pliku SVG
        int off = 50;           // Marginesy od krawędzi obrazka

        // Zapis nagłówków specyfikacji SVG
        file << "<svg width=\"" << w << "\" height=\"" << h << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"#fdfdfd\" />\n";

        // 1. RYSOWANIE ŁĄCZY RADIOWYCH (Network Backbone)
        // Rysujemy je pod spodem, jako szare przerywane linie
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double val = 0;
                try { for (int k = 0; k < K; k++) val += cplex.getValue(r[k][i * N + j]); }
                catch (...) {}

                if (val > 0.5) { // Jeśli łącze istnieje
                    file << "<line x1=\"" << (off + pos[i][0] * s) << "\" y1=\"" << (h - off - pos[i][1] * s)
                        << "\" x2=\"" << (off + pos[j][0] * s) << "\" y2=\"" << (h - off - pos[j][1] * s)
                        << "\" stroke=\"#b3b3b3\" stroke-width=\"3\" stroke-dasharray=\"6,6\"/>\n";
                }
            }
        }

        // 2. RYSOWANIE TRAJEKTORII LOTU DRONÓW
        std::string cols[] = { "red", "green", "blue", "orange", "purple", "brown" ,"yellow","pink" };
        for (int k = 0; k < K; k++) {

            // MECHANIZM PRZESUNIĘCIA WEKTOROWEGO
            // Jeżeli wiele dronów leci z punktu A do B po tej samej linii prostej, zlałyby się w jeden wektor.
            // Aby pokazać lot w formacji, przesuwamy graficznie każdego drona prostopadle o kilka pikseli.
            double dx = (k % 2 == 0 ? 1 : -1) * (k / 2 + 1) * 3.0;
            double dy = (k % 2 == 0 ? -1 : 1) * (k / 2 + 1) * 3.0;

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    double v = 0;
                    try { v = cplex.getValue(x[k][i * N + j]); }
                    catch (...) {}

                    if (i != j && v > 0.5) {
                        // Obliczenie wektorów SVG dla krawędzi przelotu
                        double x1 = off + pos[i][0] * s + dx;
                        double y1 = h - off - pos[i][1] * s + dy;
                        double x2 = off + pos[j][0] * s + dx;
                        double y2 = h - off - pos[j][1] * s + dy;

                        // Rysowanie kolorowej linii trasy
                        file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2
                            << "\" stroke=\"" << cols[k % 6] << "\" stroke-width=\"2\"/>\n";

                        // DODAWANIE TEKSTU Z DŁUGOŚCIĄ KRAWĘDZI
                        // Wyliczamy środek narysowanej linii
                        double mid_x = (x1 + x2) / 2.0;
                        double mid_y = (y1 + y2) / 2.0;

                        // Zapisanie wartości długości z macierzy odległości
                        // text-anchor="middle" idealnie wyśrodkowuje tekst.
                        file << "<text x=\"" << mid_x << "\" y=\"" << mid_y - 4
                            << "\" font-family=\"Arial\" font-size=\"10\" font-weight=\"bold\" fill=\""
                            << cols[k % 6] << "\" text-anchor=\"middle\">"
                            << std::fixed << std::setprecision(1) << dist_matrix[i][j]
                            << "</text>\n";
                    }
                }
            }

            // Rysowanie kropki docelowej - miejsca w którym dron zawisa
            for (int i = 1; i < N; i++) {
                double v = 0;
                try { v = cplex.getValue(y[k][i]); }
                catch (...) {}
                if (v > 0.5) {
                    file << "<circle cx=\"" << (off + pos[i][0] * s + dx) << "\" cy=\"" << (h - off - pos[i][1] * s + dy)
                        << "\" r=\"5\" fill=\"" << cols[k % 6] << "\"/>\n";
                }
            }
        }

        // 3. RYSOWANIE WĘZŁÓW
        for (int i = 0; i < N; i++) {
            double cx = off + pos[i][0] * s;
            double cy = h - off - pos[i][1] * s;

            // Wyróżnienie celu strategicznego (Złoty) i Bazy (Czerwony)
            std::string fill = (i == 1) ? "gold" : (i == 0 ? "red" : "black");
            double r_size = (i == 1) ? 10 : 5;

            if (i == 0) {
                // Rysuj bazę jako czerwony kwadrat
                file << "<rect x=\"" << cx - 10 << "\" y=\"" << cy - 10 << "\" width=\"20\" height=\"20\" fill=\"red\"/>";
            }
            else {
                // Rysuj cel jako kółko z czarną ramką
                file << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r_size << "\" fill=\"" << fill << "\" stroke=\"black\"/>";
            }

            // Numer punktu 
            file << "<text x=\"" << cx + 10 << "\" y=\"" << cy << "\" font-size=\"11\" font-weight=\"bold\">" << i << "</text>\n";
        }

        // Zakończenie generowania pliku
        file << "</svg>\n";
        file.close();

        std::cout << "[INFO] Pomyslnie zapisano plik wektorowy: sar_mission_map.svg" << std::endl;

        // Automatyczne wywołanie otwarcia pliku przez domyślną przeglądarkę (Windows)
        system("start sar_mission_map.svg");
    }
};

// ======================================================================================
// INTERFEJS UŻYTKOWNIKA (Funkcja pomocnicza)
// Umożliwia zmianę parametrów misji.
// ======================================================================================
template <typename T>
T getParam(std::string prompt, T defaultValue) {
    std::cout << prompt << " [Wartosc domyslna: " << defaultValue << "]: ";
    std::string line;
    std::getline(std::cin, line); // Pobranie linii wejścia

    if (line.empty()) { // Jeżeli naciśnięto tylko ENTER
        return defaultValue;
    }

    std::stringstream ss(line);
    T value;
    if (ss >> value) {
        return value;
    }
    else {
        std::cout << "Blad parsowania. Uzywam wartosci domyslnej.\n";
        return defaultValue;
    }
}

// ======================================================================================
// GŁÓWNA PĘTLA PROGRAMU (Entry Point)
// ======================================================================================
int main() {
    std::cout << "===============================================================\n";
    std::cout << " Optymalizator Liniowy MILP: Multi-UAV SAR & Relay Mission\n";
    std::cout << " Wcisnij [ENTER] aby zatwierdzic wartosci domyślne.\n";
    std::cout << "===============================================================\n\n";

    // Wartości predefiniowane (tzw. Scenariusz Badawczy nr 1)
    int def_N = 30;           // Liczba wierzchołków grafu
    int def_K = 3;            // Liczba dronów
    double def_R = 35.0;      // Skuteczny zasięg radia w terenie
    double def_Tmax = 60.0; // Budżet czasu
    double def_V = 10.0;      // Stała prędkości

    // Pobieranie konfiguracji środowiska
    int N = getParam("Rozmiar mapy (Liczba celow N)", def_N);
    int K = getParam("Rozmiar floty (Liczba dronow K)", def_K);
    double R = getParam("Zasieg radia (R w metrach)", def_R);
    double Tmax = getParam("Budzet czasu baterii (Tmax w sek)", def_Tmax);
    double V = getParam("Predkosc przelotowa maszyny (V)", def_V);

    std::cout << "\n[INFO] Rozpoczynam ladowanie srodowiska dla " << K << " dronow i " << N << " celow.\n";

    try {
        //Konstruktor planowania misji
        SARMissionPlanner planner(N, K, Tmax, R, V);

        // Definicja ułożenia celów na płaszczyźnie
        planner.generate_complex_scenario();

        // Definicja problemu fizycznego na macierz modelu CPLEX
        planner.build_model();

        // Uruchomienie solvera
        planner.solve();
    }
    catch (IloException& e) {
        // Obsługa błędów natywnych biblioteki IBM
        std::cerr << "KRYTYCZNY WYJATEK SILNIKA CPLEX: " << e << "\n";
    }
    catch (std::exception& e) {
        // Obsługa błędów pamięci i architektury C++
        std::cerr << "KRYTYCZNY WYJATEK APLIKACJI C++: " << e.what() << "\n";
    }

    std::cout << "\nWcisnij [ENTER] aby zamknac okno konsoli...";
    std::cin.get();

    return 0;
}