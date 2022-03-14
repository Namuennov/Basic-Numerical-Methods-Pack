//INTERPOLACJA LAGRANGE'A
//n = liczba argumentów, tab = wskaŸnik na tablicê z argumentami i wartoœciami funkcji, x = argument, dla którego przybli¿amy warttoœæ funkcji
double interpolationLagrange(int n, double** tab, double x);

//INTERPOLACJA NEWTONA - WERSJA REKURENCYJNA
//n = liczba argumentów, tab = wskaŸnik na tablicê z argumentami i wartoœciami funkcji,
//x = argument, dla którego przybli¿amy wartoœæ funkcji
double R_interpolationNewton(int n, double** tab, double x);

//INTERPOLACJA NEWTONA
//n = liczba argumentów, tab = wskaŸnik na tablicê z argumentami i wartoœciami funkcji,
//x = argument, dla którego przybli¿amy wartoœæ funkcji
double interpolationNewton(int n, double** tab, double x);

//METODA PROSTOK¥TÓW OBLICZANIA CA£KI OZNACZANEJ;  f = wskaŸnik na funkcjê, xp = pocz¹tek przedzia³u, xk = koniec przedzia³u, n = liczba prostok¹tów
double rectanglesMethod(double(*f)(double), double xp, double xk, int n);

//METODA TRAPEZÓW OBLICZANIA CA£KI OZNACZONEJ;  f = wskaŸnik na funkcjê, xp = pocz¹tek przedzia³u, xk = koniec przedzia³u, n = liczba trapezów
double trapezoidsMethod(double(*f)(double), double xp, double xk, int n);

//Z£O¯ONA METODA SIMPSONA;  f = wskaŸnik na funkcjê, xp = pocz¹tek przedzia³u, xk = koniec przedzia³u, n = liczba parabol
double SimpsonMethod(double(*f)(double), double xp, double xk, int n);

//METHODA MONTE CARLO OBLICZANIA CA£KI OZNACZONEJ;  f = wskaŸnik na funkcjê, xp = pocz¹tek przedzia³u, xk = koniec przedzia³u, n = liczba losowanych punktów
double MonteCarloMethod(double(*f)(double), double xp, double xk, int n);

//METODA GAUSSA ROZWI¥ZYWANIA UK£ADÓW RÓWNAÑ LINIOWYCH;  A - macierz wspó³czynników równañ liniowych, b - wektor wyrazów wolnych, n - rozmiar
double* GaussElimination(double** A, double* b, int n);

//METODA BISEKCJI ROZWI¥ZYWANIA RÓWNAÑ NIELINIOWYCH; f - wskaŸnik na funkcjê, dla której liczone jest przybli¿enie,
//xp - pocz¹tek przedzia³u, xk - koniec przedzia³u, e - dok³adnoœæ
double bisection(double(*f)(double), double xp, double xk, double e);

//REKURENCYJNA WERSJA METODY NEWTONA-RAPHSONA --> BRAK WALIDACJI DANYCH WEJŒCIOWYCH!!!
double R_NewtonRaphsonMethod(double(*f)(double), double(*f_)(double), double xp, double xk, double epsylun);

//METODA NEWTONA-RAPHSONA; f - wskaŸnik na funkcjê, dla której liczone jest przybli¿enie, f_ - wskaŸnik na pochodn¹ funkcji f,
//xp - pocz¹tek przedzia³u, xk - koniec przedzia³u, e - dok³adnoœæ
double NewtonRaphsonMethod(double(*f)(double), double(*f_)(double), double xp, double xk, double e);

//OPTYMALIZACJA METOD¥ Z£OTEGO PODZIA£U; f - wskaŸnik na funkcjê, dla której szukane jest minimum,
//xp - pocz¹tek przedzia³u, xk - koniec przedzia³u, e - dok³adnoœæ
double optimizationGold(double(*f)(double), double xp, double xk, double e);

//APROKSYMACJA LINIOWA; argumenty: tab - tablica punktów, n - rozmiar tab;
//funkcja zwraca tablicê 3 wartoœci: result[0] - wyraz wolny, result[1] - wsp. kierunkowy, result[2] - wsp. korelacji
double* approximationLinear(double** tab, int n);

//ROZWI¥ZYWANIE RÓWNAÑ RÓ¯NICZKOWYCH METOD¥ EULERA
//argumenty funkcji: f - wskaŸnik na funkcjê f z równania y' = f(x,y(x)) ,
//a - pocz¹tek przedzia³u oraz argument z zagadnienia pocz¹tkowego , b - koniec przedzia³u ,
//y0 - wartoœæ z zagadnienia pocz¹tkowego , n - liczba kroków
double** differential_equations_Euler(double(*f)(double, double), double a, double b, double y0, int n);

//ROZWI¥ZYWANIE RÓWNAÑ RÓ¯NICZKOWYCH METOD¥ HEUNA (czyli metod¹ RUNGEGO-KUTTY II RZÊDU)
//argumenty funkcji: f - wskaŸnik na funkcjê f z równania y' = f(x,y(x)) ,
//a - pocz¹tek przedzia³u oraz argument z zagadnienia pocz¹tkowego , b - koniec przedzia³u ,
//y0 - wartoœæ z zagadnienia pocz¹tkowego , n - liczba kroków
double** differential_equations_RK2(double(*f)(double, double), double a, double b, double y0, int n);

//ROZWI¥ZYWANIE RÓWNAÑ RÓ¯NICZKOWYCH METOD¥ RUNGEGO-KUTTY IV RZÊDU
//argumenty funkcji: f - wskaŸnik na funkcjê f z równania y' = f(x,y(x)) ,
//a - pocz¹tek przedzia³u oraz argument z zagadnienia pocz¹tkowego , b - koniec przedzia³u ,
//y0 - wartoœæ z zagadnienia pocz¹tkowego , n - liczba kroków
double** differential_equations_RK4(double(*f)(double, double), double a, double b, double y0, int n);

//EIGENVECTOR POWER METHOD
//argumenty: A - wskaŸnik na macierz, n - liczba wierszy macierzy A,
//m - liczba liczba kolumn macierzy A, steps - maksymalna liczba kroków
double* eigenvectorPM(double** A, int n, int m, int steps);