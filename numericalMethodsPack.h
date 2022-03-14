//INTERPOLACJA LAGRANGE'A
//n = liczba argument�w, tab = wska�nik na tablic� z argumentami i warto�ciami funkcji, x = argument, dla kt�rego przybli�amy wartto�� funkcji
double interpolationLagrange(int n, double** tab, double x);

//INTERPOLACJA NEWTONA - WERSJA REKURENCYJNA
//n = liczba argument�w, tab = wska�nik na tablic� z argumentami i warto�ciami funkcji,
//x = argument, dla kt�rego przybli�amy warto�� funkcji
double R_interpolationNewton(int n, double** tab, double x);

//INTERPOLACJA NEWTONA
//n = liczba argument�w, tab = wska�nik na tablic� z argumentami i warto�ciami funkcji,
//x = argument, dla kt�rego przybli�amy warto�� funkcji
double interpolationNewton(int n, double** tab, double x);

//METODA PROSTOK�T�W OBLICZANIA CA�KI OZNACZANEJ;  f = wska�nik na funkcj�, xp = pocz�tek przedzia�u, xk = koniec przedzia�u, n = liczba prostok�t�w
double rectanglesMethod(double(*f)(double), double xp, double xk, int n);

//METODA TRAPEZ�W OBLICZANIA CA�KI OZNACZONEJ;  f = wska�nik na funkcj�, xp = pocz�tek przedzia�u, xk = koniec przedzia�u, n = liczba trapez�w
double trapezoidsMethod(double(*f)(double), double xp, double xk, int n);

//Z�O�ONA METODA SIMPSONA;  f = wska�nik na funkcj�, xp = pocz�tek przedzia�u, xk = koniec przedzia�u, n = liczba parabol
double SimpsonMethod(double(*f)(double), double xp, double xk, int n);

//METHODA MONTE CARLO OBLICZANIA CA�KI OZNACZONEJ;  f = wska�nik na funkcj�, xp = pocz�tek przedzia�u, xk = koniec przedzia�u, n = liczba losowanych punkt�w
double MonteCarloMethod(double(*f)(double), double xp, double xk, int n);

//METODA GAUSSA ROZWI�ZYWANIA UK�AD�W R�WNA� LINIOWYCH;  A - macierz wsp�czynnik�w r�wna� liniowych, b - wektor wyraz�w wolnych, n - rozmiar
double* GaussElimination(double** A, double* b, int n);

//METODA BISEKCJI ROZWI�ZYWANIA R�WNA� NIELINIOWYCH; f - wska�nik na funkcj�, dla kt�rej liczone jest przybli�enie,
//xp - pocz�tek przedzia�u, xk - koniec przedzia�u, e - dok�adno��
double bisection(double(*f)(double), double xp, double xk, double e);

//REKURENCYJNA WERSJA METODY NEWTONA-RAPHSONA --> BRAK WALIDACJI DANYCH WEJ�CIOWYCH!!!
double R_NewtonRaphsonMethod(double(*f)(double), double(*f_)(double), double xp, double xk, double epsylun);

//METODA NEWTONA-RAPHSONA; f - wska�nik na funkcj�, dla kt�rej liczone jest przybli�enie, f_ - wska�nik na pochodn� funkcji f,
//xp - pocz�tek przedzia�u, xk - koniec przedzia�u, e - dok�adno��
double NewtonRaphsonMethod(double(*f)(double), double(*f_)(double), double xp, double xk, double e);

//OPTYMALIZACJA METOD� Z�OTEGO PODZIA�U; f - wska�nik na funkcj�, dla kt�rej szukane jest minimum,
//xp - pocz�tek przedzia�u, xk - koniec przedzia�u, e - dok�adno��
double optimizationGold(double(*f)(double), double xp, double xk, double e);

//APROKSYMACJA LINIOWA; argumenty: tab - tablica punkt�w, n - rozmiar tab;
//funkcja zwraca tablic� 3 warto�ci: result[0] - wyraz wolny, result[1] - wsp. kierunkowy, result[2] - wsp. korelacji
double* approximationLinear(double** tab, int n);

//ROZWI�ZYWANIE R�WNA� RӯNICZKOWYCH METOD� EULERA
//argumenty funkcji: f - wska�nik na funkcj� f z r�wnania y' = f(x,y(x)) ,
//a - pocz�tek przedzia�u oraz argument z zagadnienia pocz�tkowego , b - koniec przedzia�u ,
//y0 - warto�� z zagadnienia pocz�tkowego , n - liczba krok�w
double** differential_equations_Euler(double(*f)(double, double), double a, double b, double y0, int n);

//ROZWI�ZYWANIE R�WNA� RӯNICZKOWYCH METOD� HEUNA (czyli metod� RUNGEGO-KUTTY II RZ�DU)
//argumenty funkcji: f - wska�nik na funkcj� f z r�wnania y' = f(x,y(x)) ,
//a - pocz�tek przedzia�u oraz argument z zagadnienia pocz�tkowego , b - koniec przedzia�u ,
//y0 - warto�� z zagadnienia pocz�tkowego , n - liczba krok�w
double** differential_equations_RK2(double(*f)(double, double), double a, double b, double y0, int n);

//ROZWI�ZYWANIE R�WNA� RӯNICZKOWYCH METOD� RUNGEGO-KUTTY IV RZ�DU
//argumenty funkcji: f - wska�nik na funkcj� f z r�wnania y' = f(x,y(x)) ,
//a - pocz�tek przedzia�u oraz argument z zagadnienia pocz�tkowego , b - koniec przedzia�u ,
//y0 - warto�� z zagadnienia pocz�tkowego , n - liczba krok�w
double** differential_equations_RK4(double(*f)(double, double), double a, double b, double y0, int n);

//EIGENVECTOR POWER METHOD
//argumenty: A - wska�nik na macierz, n - liczba wierszy macierzy A,
//m - liczba liczba kolumn macierzy A, steps - maksymalna liczba krok�w
double* eigenvectorPM(double** A, int n, int m, int steps);