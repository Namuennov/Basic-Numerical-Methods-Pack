#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "numericalMethods.h"

#define _USE_MATH_DEFINES
#define FI 0.61803398


template <class T> void swap(T& a, T& b)
{
    T c(a);
    a = b;
    b = c;
}

template <class T> void heapifyInterpolation(T** tab, int size, int i)
{
    for (int l = 2 * i; l <= size; i = l, l = 2 * i)
    {
        if (l + 1 <= size && tab[l + 1][0] > tab[l][0]) ++l;
        if (tab[i][0] >= tab[l][0]) return;
        swap(tab[i], tab[l]);
    }
}

template <class T> void build_heapInterpolation(T** tab, int size)
{
    for (long i = size / 2; i; --i)
        heapifyInterpolation(tab, size, i);
}

template <class T> void heap_sortInterpolation(T** tab, int size)
{
    --tab;
    build_heapInterpolation(tab, size);
    while (size > 1)
    {
        swap(tab[1], tab[size]);
        heapifyInterpolation(tab, --size, 1);
    }
}

//INTERPOLACJA LAGRANGE'A
//n = liczba argumentów, tab = wskaźnik na tablicę z argumentami i wartościami funkcji, x = argument, dla którego przybliżamy warttość funkcji
double interpolationLagrange(int n, double** tab, double x)
{
    //sprawdzenie czy rozmiar tablicy jest liczbą dodatnią
    if (n <= 0)
    {
        std::cout << "Wrong arguments! (size of table)\n";
        return 0.;
    }
    //sprawdzenie czy nie wysłano pustego wskaźnika
    if (tab == nullptr)
    {
        std::cout << "Wrong arguments! (empty pointer to table)\n";
        return 0.;
    }

    //przepisanie tablicy wejściowej do tablicy Stab - jest to potrzebne, żeby nie sortować tablicy wejściowej
    double** Stab = new double* [n];
    for (int i = 0; i < n; i++) Stab[i] = new double[2];
    for (int i = 0; i < n; i++)
    {
        Stab[i][0] = tab[i][0];
        Stab[i][1] = tab[i][1];
    }
    //posortowanie tablicy Stab
    heap_sortInterpolation(Stab, n);

    //sprawdzenie czy x jest w zadanym przedziale
    if (x < Stab[0][0] || x > Stab[n - 1][0])
    {
        std::cout << "Wrong arguments!\nArgument, for which you want to approximate value, is out of numerical interval sent in table\n";
        return 0.;
    }

    //sprawdzenie czy każdy x w tablicy wejścowej jest unikalny
    //poniższa pętla przerzuca komórki powtarzających się x'ów na początek tablicy oraz ustawia parametry
    //do liczenia (będących niżej) tablicy l[] i zmiennej result, czyli ,,wskaźnik" na początek tablicy, dla
    //której liczymy (zmienna ptr) i rozmiar tej tablicy (zmodyfikowana zmienna wejściowa n)
    int m = n - 1;
    int ptr = 0;
    for (int i = 0; i < m; i++)
    {
        if (Stab[i][0] == Stab[i + 1][0])
        {
            swap(Stab[i], Stab[ptr]);
            ptr++;
            n--;
        }
    }
    //przesunięcie początku tablicy o wyliczoną wartość
    Stab += ptr;



    //zaalokowanie pamięci potrzebnej na tablicę ze współczynnikami l
    double* l = new double[n];
    //wyznaczenie współczynników l
    for (int i = 0; i < n; i++)
    {
        //zmienna z aktualnym współczynnikiem - ponieważ będzie na niej wykonywana operacja *=, to jej początkowa wartość =1
        double l_i = 1.;
        //wyznaczenie aktualnego współczynnika w wewnętrznej pętli
        for (int j = 0; j < n; j++)
        {
            //podstawienie do wzoru
            if (j != i) l_i *= ((x - Stab[j][0])) / ((Stab[i][0]) - (Stab[j][0]));
            //jeżeli dotarliśmy do wierzchołka, dla którego liczymy l, to musimy pominąć iterację
            else continue;
        }
        //zapisanie wartości zmiennej do tablicy
        l[i] = l_i;
    }

    //utworzenie zmiennej na wynik
    double result = 0.;
    //obliczenie wyniku przy pomocy wyliczonych wcześniej współczynników l
    for (int i = 0; i < n; i++)
    {
        //podstawienie do wzoru
        result += Stab[i][1] * l[i];
    }



    //przesunięcie wskaźnika Stab na poprzednią wartość
    Stab -= ptr;
    //dealokacja pamięci
    delete[] l;
    for (int i = 0; i < n; i++) delete[] * (Stab + i);
    delete[] Stab;
    //zwrócenie wyniku
    return result;
}

//rekurencyjna funkcja do obliczania ilorazów różnicowych
//n = liczba argumentów, tab = wskaźnik na tablicę z argumentami i wartościami funkcji
double newtonFactor(int n, double** tab)
{
    //walidacja danych - wskaźnik
    if (tab == nullptr)
    {
        std::cout << "Wrong arguments!\n";
        return 0.;
    }

    //rekurencyjne obliczenie ilorazów różnicowych
    if (n > 1) return ((newtonFactor((n - 1), (tab + 1)) - newtonFactor((n - 1), tab)) / (tab[n - 1][0] - tab[0][0]));
    else if (n == 1) return tab[0][1];
    //walidacja danych - rozmiar tablicy
    else
    {
        std::cout << "Wrong arguments!\n";
        return 0.;
    }
}

//INTERPOLACJA NEWTONA - WERSJA REKURENCYJNA
//n = liczba argumentów, tab = wskaźnik na tablicę z argumentami i wartościami funkcji,
//x = argument, dla którego przybliżamy wartość funkcji
double R_interpolationNewton(int n, double** tab, double x)
{
    //sprawdzenie czy rozmiar tablicy jest liczbą dodatnią
    if (n <= 0)
    {
        std::cout << "Wrong arguments! (size of table)\n";
        return 0.;
    }
    //sprawdzenie czy nie wysłano pustego wskaźnika
    if (tab == nullptr)
    {
        std::cout << "Wrong arguments! (empty pointer to table)\n";
        return 0.;
    }

    //przepisanie tablicy wejściowej do tablicy Stab - jest to potrzebne, żeby nie sortować tablicy wejściowej
    double** Stab = new double* [n];
    for (int i = 0; i < n; i++) Stab[i] = new double[2];
    for (int i = 0; i < n; i++)
    {
        Stab[i][0] = tab[i][0];
        Stab[i][1] = tab[i][1];
    }
    //posortowanie tablicy Stab
    heap_sortInterpolation(Stab, n);

    //sprawdzenie czy x jest w zadanym przedziale
    if (x < Stab[0][0] || x > Stab[n - 1][0])
    {
        std::cout << "Wrong arguments!\nArgument, for which you want to approximate value, is out of numerical interval sent in table\n";
        return 0.;
    }

    //sprawdzenie czy każdy x w tablicy wejścowej jest unikalny
    //poniższa pętla przerzuca komórki powtarzających się x'ów na początek tablicy oraz ustawia parametry
    //do liczenia (będących niżej) tablicy l[] i zmiennej result, czyli ,,wskaźnik" na początek tablicy, dla
    //której liczymy (zmienna ptr) i rozmiar tej tablicy (zmodyfikowana zmienna wejściowa n)

    //pomocnicza zmienna m
    int m = n - 1;
    int ptr = 0;
    for (int i = 0; i < m; i++)
    {
        if (Stab[i][0] == Stab[i + 1][0])
        {
            swap(Stab[i], Stab[ptr]);
            ptr++;
            n--;
        }
    }
    //przesunięcie początku tablicy o wyliczoną wartość
    Stab += ptr;

    //zaalokowanie tablicy na ilorazy różnicowe
    double* l = new double[n];
    //obliczenie ilorazów różnicowych funkcją rekurencyjną - pętla wpisuje do tablicy tylko ilorazy wykorzystywane do obliczenia przybliżenia
    for (int i = 0; i < n; i++) l[i] = newtonFactor((i + 1), Stab);

    //zaalokowanie tablicy na ,,współczynniki" od (x-x0) do (x-x0)(x-x1)...(x-xn-1)
    double* l2 = new double[n];
    //sztuczny współczynnik równy 1, żeby w następnej pętli wyszedł prawidłowy wynik
    l2[0] = 1.;
    //obliczenie ,,współczynników"
    for (int i = 1; i < n; i++) l2[i] = (l2[i - 1] * (x - Stab[i - 1][0]));

    //obliczenie przybliżenia dla x'a zgodnie z wzorem
    double result = 0.;
    for (int i = 0; i < n; i++) result += (l[i] * l2[i]);

    //przesunięcie wskaźnika Stab na poprzednią wartość
    Stab -= ptr;
    //dealokacja pamięci
    delete[] l;
    delete[] l2;
    for (int i = 0; i < n; i++) delete[] Stab[i];
    delete[] Stab;
    //zwrócenie wyniku
    return result;
}

//INTERPOLACJA NEWTONA
//n = liczba argumentów, tab = wskaźnik na tablicę z argumentami i wartościami funkcji,
//x = argument, dla którego przybliżamy wartość funkcji
double interpolationNewton(int n, double** tab, double x)
{
    //sprawdzenie czy rozmiar tablicy jest liczbą dodatnią
    if (n <= 0)
    {
        std::cout << "Wrong arguments! (size of table)\n";
        return 0.;
    }
    //sprawdzenie czy nie wysłano pustego wskaźnika
    if (tab == nullptr)
    {
        std::cout << "Wrong arguments! (empty pointer to table)\n";
        return 0.;
    }

    //przepisanie tablicy wejściowej do tablicy Stab - jest to potrzebne, żeby nie sortować tablicy wejściowej
    double** Stab = new double* [n];
    for (int i = 0; i < n; i++) Stab[i] = new double[2];
    for (int i = 0; i < n; i++)
    {
        Stab[i][0] = tab[i][0];
        Stab[i][1] = tab[i][1];
    }
    //posortowanie tablicy Stab
    heap_sortInterpolation(Stab, n);

    //sprawdzenie czy x jest w zadanym przedziale
    if (x < Stab[0][0] || x > Stab[n - 1][0])
    {
        std::cout << "Wrong arguments!\nArgument, for which you want to approximate value, is out of numerical interval sent in table\n";
        return 0.;
    }

    //sprawdzenie czy każdy x w tablicy wejścowej jest unikalny
    //poniższa pętla przerzuca komórki powtarzających się x'ów na początek tablicy oraz ustawia parametry
    //do liczenia (będących niżej) tablicy l[] i zmiennej result, czyli ,,wskaźnik" na początek tablicy, dla
    //której liczymy (zmienna ptr) i rozmiar tej tablicy (zmodyfikowana zmienna wejściowa n)

    //pomocnicza zmienna m
    int m = n - 1;
    int ptr = 0;
    for (int i = 0; i < m; i++)
    {
        if (Stab[i][0] == Stab[i + 1][0])
        {
            swap(Stab[i], Stab[ptr]);
            ptr++;
            n--;
        }
    }
    //przesunięcie początku tablicy o wyliczoną wartość
    Stab += ptr;

    //zmienna pomocnicza m przyspieszy obliczenia, ponieważ n-1 występuje w funkcji wiele razy - trzeba zaktualizować jej wartość
    m = n - 1;
    //zaalokowanie tablicy na ilorazy różnicowe
    double** l = new double* [m];
    //przesunięcie wskaźnika w celu uniknięcia zbędnego dodawanie jedynki
    l--;
    //zaalokowanie drugiego wymiaru tablicy na ilorazy różnicowe
    for (int i = m; i > 0; i--) l[n - i] = new double[i];
    //przesunięcie wskaźnika spowrotem
    l++;

    //obliczenie ilorazów różnicowych ,,w pierwszej kolumnie" - oblicza się je inaczej niż te następne, bo podstawia się y, a nie inne ilorazy
    for (int i = 0; i < m; i++) l[0][i] = ((Stab[i + 1][1] - Stab[i][1]) / (Stab[i + 1][0] - Stab[i][0]));
    //obliczenie reszty ilorazów różnicowych
    for (int i = 1; i < m; i++) for (int j = 0; j < (m - i); j++) l[i][j] = (l[i - 1][j + 1] - l[i - 1][j]) / (Stab[j + i + 1][0] - Stab[j][0]);

    //zaalokowanie tablicy na ,,współczynniki" od (x-x0) do (x-x0)(x-x1)...(x-xn-1)
    double* l2 = new double[n];
    //sztuczny współczynnik równy 1, żeby w następnej pętli wyszedł prawidłowy wynik
    l2[0] = 1.;
    //obliczenie ,,współczynników"
    for (int i = 1; i < n; i++) l2[i] = (l2[i - 1] * (x - Stab[i - 1][0]));

    //obliczenie przybliżenia dla x'a zgodnie z wzorem
    double result = Stab[0][1];
    for (int i = 0; i < m; i++) result += (l[i][0] * l2[i + 1]);

    //przesunięcie wskaźnika Stab na poprzednią wartość
    Stab -= ptr;
    //dealokacja pamięci
    for (int i = 0; i < m; i++) delete[] l[i];
    delete[] l;
    delete[] l2;
    for (int i = 0; i < n; i++) delete[] Stab[i];
    delete[] Stab;
    //zwrócenie wyniku
    return result;
}

//METODA PROSTOKĄTÓW OBLICZANIA CAŁKI OZNACZANEJ;  f = wskaźnik na funkcję, xp = początek przedziału, xk = koniec przedziału, n = liczba prostokątów
double rectanglesMethod(double(*f)(double), double xp, double xk, int n)
{
    //walidacja argumentów
    if (xk <= xp)
    {
        std::cout << "Wrong arguments! (wrong numerical interval)\n";
        return 0.;
    }
    if (n <= 0)
    {
        std::cout << "Wrong arguments! (amount of rectangles)\n";
        return 0.;
    }
    if (f == nullptr)
    {
        std::cout << "Wrong arguments! (pointer to function)\n";
        return 0.;
    }

    //obliczenie dx
    double dx = (xk - xp) / n;
    //,,selekcja" argumentów
    //alokacja pamięci
    double* x_tab = new double[n];
    //znalezienie potrzebnych x'ów
    x_tab[0] = xp + dx;
    for (int i = 1; i < (n - 1); i++) x_tab[i] = x_tab[i - 1] + dx;
    //lepiej x_tab[n-1] przypisać tak, bo przez niedoskonałość obliczeń komputerowych może wyjść nieco inna liczba
    x_tab[n - 1] = xk;
    //obliczenie wyniku
    double result = 0.;
    for (int i = 0; i < n; i++) result += f(x_tab[i]);
    //można ,,wyciągnąć dx przed nawias" i pomnożyć tylko raz
    result *= dx;
    //dealokacja pamięci
    delete[] x_tab;
    //zwrócenie wyniku
    return result;
}

//METODA TRAPEZÓW OBLICZANIA CAŁKI OZNACZONEJ;  f = wskaźnik na funkcję, xp = początek przedziału, xk = koniec przedziału, n = liczba trapezów
double trapezoidsMethod(double(*f)(double), double xp, double xk, int n)
{
    //walidacja argumentów
    if (xk <= xp)
    {
        std::cout << "Wrong arguments! (wrong numerical interval)\n";
        return 0.;
    }
    if (n <= 0)
    {
        std::cout << "Wrong arguments! (amount of trapezoids)\n";
        return 0.;
    }
    if (f == nullptr)
    {
        std::cout << "Wrong arguments! (pointer to function)\n";
        return 0.;
    }

    //obliczenie dx
    double dx = (xk - xp) / n;
    //,,selekcja" argumentów
    //alokacja pamięci
    double* x_tab = new double[n + 1];
    //znalezienie potrzebnych x'ów
    x_tab[0] = xp;
    for (int i = 1; i < n; i++) x_tab[i] = x_tab[i - 1] + dx;
    //lepiej x_tab[n-1] przypisać tak, bo przez niedoskonałość obliczeń komputerowych może wyjść nieco inna liczba
    x_tab[n] = xk;
    //tablica y'ów, żeby nie powielać obliczeń
    double* y_tab = new double[n + 1];
    for (int i = 0; i < (n + 1); i++) y_tab[i] = f(x_tab[i]);
    //obliczenie wyniku
    //zsumowanie składników sumy, które pojawiają się tylko raz
    double result = y_tab[0] + y_tab[n];
    //dodanie do sumy składników pojawiających się 2 razy
    for (int i = 1; i < n; i++) result += (y_tab[i] * 2);
    //można ,,wyciągnąć dx/2 przed nawias" i pomnożyć tylko raz
    result *= (dx / 2);
    //dealokacja pamięci
    delete[] x_tab;
    delete[] y_tab;
    //zwrócenie wyniku
    return result;
}

//funkcja pomocnicza do losowania pseudolosowych liczb typu double z przedziału <a;b)
double d_rand(double a, double b)
{
    double d = (double)rand() / RAND_MAX;
    return (a + d * (b - a));
}

//ZŁOŻONA METODA SIMPSONA;  f = wskaźnik na funkcję, xp = początek przedziału, xk = koniec przedziału, n = liczba parabol
double SimpsonMethod(double(*f)(double), double xp, double xk, int n)
{
    //walidacja argumentów
    if (xk <= xp)
    {
        std::cout << "Wrong arguments! (wrong numerical interval)\n";
        return 0.;
    }
    if (n <= 0)
    {
        std::cout << "Wrong arguments! (amount of parabolas)\n";
        return 0.;
    }
    if (f == nullptr)
    {
        std::cout << "Wrong arguments! (pointer to function)\n";
        return 0.;
    }

    //utworzenie zmiennej pomocniczej m, ponieważ n*2 często się pojawia
    int m = n * 2;
    //obliczenie h
    double h = (xk - xp) / m;
    //utworzenie tablicy dynamicznej z wartościami funkcji, dzięki czemu nie trzeba powielać wielu obliczeń
    double* Ytab = new double[m + 1];
    for (int i = 0; i < m; i++) Ytab[i] = f((xp + (h * i)));
    //w celu zwiększenia dokładności ostatni element tablicy liczony jest poza pętlą
    Ytab[m] = f(xk);
    //obliczenie wyniku
    double result = 0.;
    for (int i = 0; i < n; i++)
    {
        int j = (i * 2);
        result += (Ytab[j] + (4 * Ytab[(j + 1)]) + Ytab[(j + 2)]);
    }
    result *= (h / 3.);
    //dealokacja pamięci
    delete[] Ytab;
    return result;
}

//METHODA MONTE CARLO OBLICZANIA CAŁKI OZNACZONEJ;  f = wskaźnik na funkcję, xp = początek przedziału, xk = koniec przedziału, n = liczba losowanych punktów
double MonteCarloMethod(double(*f)(double), double xp, double xk, int n)
{
    //walidacja argumentów
    if (xk <= xp)
    {
        std::cout << "Wrong arguments! (wrong numerical interval)\n";
        return 0.;
    }
    if (n <= 0)
    {
        std::cout << "Wrong arguments! (amount of drawn points)\n";
        return 0.;
    }
    if (f == nullptr)
    {
        std::cout << "Wrong arguments! (pointer to function)\n";
        return 0.;
    }

    //obliczenie wartości funkcji dla losowych elementów z przedziału przy użyciu pomocniczej funkcji d_rand()
    double Favg = 0.;
    for (int i = 0; i < n; i++) Favg += f(d_rand(xp, xk));
    //obliczenie xk-xp i wartości bezwzględnej z tego - dzięki takiemu liczeniu nie trzeba załączać dodatkowej biblioteki z abs()
    double X = xk - xp;
    if (X < 0.) X *= (-1);
    //obliczenie i zwrócenie wyniku
    Favg *= (X / n);
    return Favg;
}

//METODA GAUSSA ROZWIĄZYWANIA UKŁADÓW RÓWNAŃ LINIOWYCH;  A - macierz współczynników równań liniowych, b - wektor wyrazów wolnych, n - rozmiar
double* GaussElimination(double** A, double* b, int n)
{
    //walidacja danych wejściowych
    if (A == nullptr)
    {
        printf("Wrong arguments! (matrix A)\n");
        return nullptr;
    }
    if (b == nullptr)
    {
        printf("Wrong arguments! (vector b)\n");
        return nullptr;
    }
    if (n <= 0)
    {
        printf("Wrong arguments! (size lesser than 0)\n");
        return nullptr;
    }
    if (A[n - 1][n - 1] == NULL)
    {
        printf("Wrong arguments! (size of matrix A)\n");
        return nullptr;
    }
    if (b[n - 1] == NULL)
    {
        printf("Wrong arguments! (size of vector b)\n");
        return nullptr;
    }
    if (A[0][0] == 0.)
    {
        printf("Wrong arguments! (first element of matrix A equals 0)\n");
        return nullptr;
    }

    //alokacja pamięci na tablicę pomocniczą temp[][] - jako że ,,użytkownik nie widzi",
    //co się dzieje w funkcji, to można wpisać macierz A i wektor b do jednej tablicy
    double** temp = new double* [n];
    for (int i = 0; i < n; i++) temp[i] = new double[n + 1];

    //alokacja tablicy na wynik
    double* result = new double[n];

    //przepisanie do utworzonej tablicy wartości z tablic wejściowych
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) temp[i][j] = A[i][j];
    for (int i = 0; i < n; i++) temp[i][n] = b[i];

    //wyliczenie wartości znajdujących się w macierzy trójkątnej górnej
    //pierwsza pętla (ze zmienną i) odpowiada ,,etapowi", czyli eliminacji zmiennej Xi z równań
    for (int i = 0; i < (n - 1); i++)
    {
        //liczba przez którą dzielimy - wspólna dla całego etapu
        double a = temp[i][i];
        //druga pętla (ze zmienną j) odpowiada wierszowi, który jest modyfikowany
        for (int j = (n - 1); j > i; j--)
        {
            //współczynnik potrzebny do liczenia wartości elementu - wspólny dla całego wiersza
            double b = (temp[j][i] / a);
            //trzecia pętla (ze zmienną k) odpowiada elementom wiersza modyfikowanego
            for (int k = (i + 1); k < (n + 1); k++) temp[j][k] -= (b * temp[i][k]);
        }
    }

    //obliczenie wartości niewiadomych, czyli elementów tabeli result[]
    //pierwsza pętla (ze zmienną i) odpowiada liczonej niewiadomej Xi
    for (int i = (n - 1); i >= 0; i--)
    {
        //przepisanie wartości ze ,,zmodyfikowanego wektora b"
        result[i] = temp[i][n];
        //druga pętla odpowiada za odejmowanie kolejnych wartości od wyrazu wolnego
        for (int j = (n - i - 1); j > 0; j--) result[i] -= (temp[i][j + i] * result[j + i]);
        //podzielenie przez liczbę w celu uzyskania ostatecznego wyniku
        result[i] /= temp[i][i];
    }

    //dealokacja pamięci
    for (int i = 0; i < n; i++) delete[] temp[i];
    delete[] temp;
    //zwrócenie wyniku
    return result;
}

//METODA BISEKCJI ROZWIĄZYWANIA RÓWNAŃ NIELINIOWYCH; f - wskaźnik na funkcję, dla której liczone jest przybliżenie,
//xp - początek przedziału, xk - koniec przedziału, e - dokładność
double bisection(double(*f)(double), double xp, double xk, double e)
{
    //walidacja danych
    if (xk <= xp)
    {
        printf("Wrong arguments! (wrong numerical interval)\n");
        return NULL;
    }
    if (f == nullptr)
    {
        printf("Wrong arguments! (pointer to function)\n");
        return NULL;
    }
    //sprawdzenie warunku f(xk)*f(xp) < 0
    if ((f(xk) * f(xp)) >= 0.)
    {
        if ((f(xk) * f(xp)) > 0.)
        {
            printf("Wrong arguments! (f(xp)*f(xk) > 0)\n");
            return NULL;
        }
        //jeśli początek albo koniec przedziału jest argumentem zerowym, to zwrócenie tego argumentu
        else if (f(xp) == 0.) return xp;
        else return xk;
    }

    //zapisanie do zmiennej result środkowego argumentu
    double result = (xk + xp) / 2;
    //zmienna pomocnicza inicjalizowana tutaj, żeby nie trzeba było jej inicjalizować co iterację pętli
    double d;
    //pętla nieskończona - zakończy się dopiero po zwróceniu wyniku
    while (1)
    {
        //sprawdzenie czy (w wypadku nie znalezienia wcześniej dokładnego pierwiastka) nie osiągnięto zadanej dokładności
        if (abs(xk - xp) <= e) return result;
        //do zmiennej pomocniczej d zapisywana jest wartość f(result), ponieważ wykorzystuje się ją kilkukrotnie
        d = f(result);
        //jeśli znaleziono pierwiastek, to zwrócenie wartości
        if (d == 0) return result;
        //sprawdzenie, dla której połowy przedziału trzeba teraz liczyć
        else if (d * f(xp) < 0.) xk = result;
        else if (d * f(xk) < 0.) xp = result;
        //zapisanie do zmiennej result środkowego argumentu
        result = (xk + xp) / 2;
    }
}

//funkcja pomocnicza do rekurencyjnej wersji Metody Newtona-Raphsona
double NR(double(*f)(double), double(*f_)(double), double xn, double epsylun)
{
    double a = f(xn);
    if (abs(a) <= epsylun) return xn;
    double xn_ = (xn - (a / f_(xn)));
    return NR(f, f_, xn_, epsylun);
}

//REKURENCYJNA WERSJA METODY NEWTONA-RAPHSONA --> BRAK WALIDACJI DANYCH WEJŚCIOWYCH!!!
double R_NewtonRaphsonMethod(double(*f)(double), double(*f_)(double), double xp, double xk, double epsylun)
{
    return NR(f, f_, xk, epsylun);
}

//METODA NEWTONA-RAPHSONA; f - wskaźnik na funkcję, dla której liczone jest przybliżenie, f_ - wskaźnik na pochodną funkcji f,
//xp - początek przedziału, xk - koniec przedziału, e - dokładność
double NewtonRaphsonMethod(double(*f)(double), double(*f_)(double), double xp, double xk, double e)
{
    //walidacja danych
    if (xk <= xp)
    {
        printf("Wrong arguments! (wrong numerical interval)\n");
        return NULL;
    }
    if (f == nullptr)
    {
        printf("Wrong arguments! (pointer to function)\n");
        return NULL;
    }
    if (f_ == nullptr)
    {
        printf("Wrong arguments! (pointer to derivative of a function)\n");
        return NULL;
    }
    //sprawdzenie warunku f(xk)*f(xp) < 0
    if ((f(xk) * f(xp)) >= 0.)
    {
        if ((f(xk) * f(xp)) > 0.)
        {
            printf("Wrong arguments! (f(xp)*f(xk) > 0)\n");
            return NULL;
        }
        //jeśli początek albo koniec przedziału jest argumentem zerowym, to zwrócenie tego argumentu
        else if (f(xp) == 0.) return xp;
        else return xk;
    }

    //zmienna pomocnicza inicjalizowana tutaj, żeby nie trzeba było jej inicjalizować co iterację pętli
    double a;
    //pętla nieskończona - zakończy się dopiero po spełnieniu kryterium stopu
    while (1)
    {
        //do zmiennej pomocniczej a zapisywana jest wartość f(xk), ponieważ zazwyczaj wykorzystuje się ją 2 razy
        a = f(xk);
        //sprawdzenie czy spełnione jest kryterium stopu - jeśli tak, to zwrócenie wyniku
        if ((abs(a) <= e) || (abs(xk - xp) <= e)) return xk;
        //xp przyjmuje wartość xn ze wzoru
        xp = xk;
        //xk przyjmuje wartość x(n+1) ze wzoru
        xk = (xk - (a / f_(xk)));
    }
}

//OPTYMALIZACJA METODĄ ZŁOTEGO PODZIAŁU; f - wskaźnik na funkcję, dla której szukane jest minimum,
//xp - początek przedziału, xk - koniec przedziału, e - dokładność
double optimizationGold(double(*f)(double), double xp, double xk, double e)
{
    //walidacja danych
    if (xk <= xp)
    {
        printf("Wrong arguments! (wrong numerical interval)\n");
        return NULL;
    }
    if (f == nullptr)
    {
        printf("Wrong arguments! (pointer to function)\n");
        return NULL;
    }
    if (e < 0)
    {
        printf("Wrong arguments! (precision is below 0)\n");
        return NULL;
    }

    //zmienne pomocnicze inicjalizowane tutaj, żeby nie trzeba było ich inicjalizować co iterację pętli
    double xl, xr, diff;
    //pętla nieskończona - zakończy się dopiero po zwróceniu wyniku
    while (1)
    {
        //obliczenie różnicy xk-xp, ponieważ jest wykorzystywana 3 razy w 1 iteracji
        diff = xk - xp;
        //zwrócenie wyniku w przypadku osiągnięcia zadanej dokładności
        if (abs(diff) <= e) return (xk - (diff / 2.));

        //obliczenie xl (x left) i xr (x right)
        xl = xk - (FI * diff);
        xr = xp + (FI * diff);
        //wyznaczenie początku/końca przedziału do następnej iteracji
        if (f(xl) > f(xr)) xp = xl;
        else xk = xr;
    }
}

//APROKSYMACJA LINIOWA; argumenty: tab - tablica punktów, n - rozmiar tab;
//funkcja zwraca tablicę 3 wartości: result[0] - wyraz wolny, result[1] - wsp. kierunkowy, result[2] - wsp. korelacji
double* approximationLinear(double** tab, int n)
{
    //walidacja danych wejściowych
    if (tab == nullptr)
    {
        printf("Bad arguments! (pointer to table)\n");
        return nullptr;
    }
    if (n < 2)
    {
        printf("Bad arguments! (size of table)\n");
        return nullptr;
    }
    //sprawdzenie czy nie x = const (wtedy 1 punkt lub pionowa linia)
    bool b = false;
    for (int i = 1; i < n; i++)
    {
        if (tab[i][0] != tab[0][0])
        {
            b = true;
            break;
        }
    }
    if (b == false)
    {
        printf("Bad parametres! (you tried to find vertical line or whole table contains 1 repeated point)\n");
        return nullptr;
    }

    //tablica pomocnicza na sumy powtarzające się we wzorach
    double sum[5] = { 0.,0.,0.,0.,0. };
    //sum[0] = suma x'ów
    for (int i = 0; i < n; i++) sum[0] += tab[i][0];
    //sum[1] = suma y'ów
    for (int i = 0; i < n; i++) sum[1] += tab[i][1];
    //sum[2] = suma iloczynów x,y
    for (int i = 0; i < n; i++) sum[2] += (tab[i][0] * tab[i][1]);
    //sum[3] = suma kwadratów x'ów
    for (int i = 0; i < n; i++) sum[3] += (tab[i][0] * tab[i][0]);
    //sum[4] = suma kwadratów y'ów
    for (int i = 0; i < n; i++) sum[4] += (tab[i][1] * tab[i][1]);
    //tablica na wynik
    double* result = new double[3];
    //podstawienie do wzorów
    result[0] = ((sum[1] * sum[3]) - (sum[0] * sum[2])) / ((n * sum[3]) - (sum[0] * sum[0]));
    result[1] = ((n * sum[2]) - (sum[0] * sum[1])) / ((n * sum[3]) - (sum[0] * sum[0]));
    //przypadek dla funkcji stałej
    if (result[1] == 0.) result[2] = 0.;
    else result[2] = ((n * sum[2]) - (sum[0] * sum[1])) / (sqrt((n * sum[3]) - (sum[0] * sum[0])) * sqrt((n * sum[4]) - (sum[1] * sum[1])));
    //zwrócenie wyniku
    return result;
}

//ROZWIĄZYWANIE RÓWNAŃ RÓŻNICZKOWYCH METODĄ EULERA
//argumenty funkcji: f - wskaźnik na funkcję f z równania y' = f(x,y(x)) ,
//a - początek przedziału oraz argument z zagadnienia początkowego , b - koniec przedziału ,
//y0 - wartość z zagadnienia początkowego , n - liczba kroków
double** differential_equations_Euler(double(*f)(double, double), double a, double b, double y0, int n)
{
    //walidacja argumentów
    if (f == nullptr)
    {
        printf("Wrong arguments! (pointer to function)\n");
        return nullptr;
    }
    if (b <= a)
    {
        printf("Wrong arguments! (numerical interval)\n");
        return nullptr;
    }
    if (n <= 0)
    {
        printf("Wrong arguments! (amount of steps)\n");
        return nullptr;
    }

    //obliczenie długości kroku
    double h = (b - a) / static_cast<double>(n);
    //zaalokowanie zmiennej na wynik
    double** result = new double* [n];
    for (int i = 0; i < n; i++) result[i] = new double[2];
    //obliczenie współrzędnych punktów
    result[0][0] = a + h;
    result[0][1] = y0 + (h * f(a, y0));
    n--;
    //zmienne pomocnicze do zmniejszenia liczby odczytów z pamięci
    double px, py;
    for (int i = 1; i < n; i++)
    {
        px = result[i - 1][0];
        py = result[i - 1][1];
        result[i][0] = px + h;
        result[i][1] = py + (h * f(px, py));
    }
    //przypisanie ,,na sztywno" wartości ostatniej odciętej
    result[n][0] = b;
    if (n >= 1) result[n][1] = result[n - 1][1] + (h * f(result[n - 1][0], result[n - 1][1]));
    //zwrócenie wyniku
    return result;
}

//ROZWIĄZYWANIE RÓWNAŃ RÓŻNICZKOWYCH METODĄ HEUNA (czyli metodą RUNGEGO-KUTTY II RZĘDU)
//argumenty funkcji: f - wskaźnik na funkcję f z równania y' = f(x,y(x)) ,
//a - początek przedziału oraz argument z zagadnienia początkowego , b - koniec przedziału ,
//y0 - wartość z zagadnienia początkowego , n - liczba kroków
double** differential_equations_RK2(double(*f)(double, double), double a, double b, double y0, int n)
{
    //walidacja argumentów
    if (f == nullptr)
    {
        printf("Wrong arguments! (pointer to function)\n");
        return nullptr;
    }
    if (b <= a)
    {
        printf("Wrong arguments! (numerical interval)\n");
        return nullptr;
    }
    if (n <= 0)
    {
        printf("Wrong arguments! (amount of steps)\n");
        return nullptr;
    }

    //obliczenie długości kroku
    double h = (b - a) / static_cast<double>(n);
    //zaalokowanie zmiennej na wynik
    double** result = new double* [n];
    for (int i = 0; i < n; i++) result[i] = new double[2];
    //zmienne pomocnicze do zwiększenia wydajności
    double px = a + h;
    double py;
    double pr = f(a, y0);
    //obliczenie współrzędnych punktów
    result[0][0] = px;
    result[0][1] = y0 + ((h / 2.) * (pr + f(a, y0 + (h * pr))));
    n--;
    for (int i = 1; i < n; i++)
    {
        py = result[i - 1][1];
        pr = f(px, py);
        result[i][1] = py + ((h / 2.) * (pr + f(px, py + (h * pr))));
        px += h;
        result[i][0] = px;
    }
    //przypisanie ,,na sztywno" wartości współrzędnych ostatniego punktu
    result[n][0] = b;
    if (n >= 1) result[n][1] = result[n - 1][1] + ((h / 2.) * (f(b - h, result[n - 1][1]) + f(b - h, result[n - 1][1] + (h * f(b - h, result[n - 1][1])))));
    //zwrócenie wyniku
    return result;
}

//ROZWIĄZYWANIE RÓWNAŃ RÓŻNICZKOWYCH METODĄ RUNGEGO-KUTTY IV RZĘDU
//argumenty funkcji: f - wskaźnik na funkcję f z równania y' = f(x,y(x)) ,
//a - początek przedziału oraz argument z zagadnienia początkowego , b - koniec przedziału ,
//y0 - wartość z zagadnienia początkowego , n - liczba kroków
double** differential_equations_RK4(double(*f)(double, double), double a, double b, double y0, int n)
{
    //walidacja argumentów
    if (f == nullptr)
    {
        printf("Wrong arguments! (pointer to function)\n");
        return nullptr;
    }
    if (b <= a)
    {
        printf("Wrong arguments! (numerical interval)\n");
        return nullptr;
    }
    if (n <= 0)
    {
        printf("Wrong arguments! (amount of steps)\n");
        return nullptr;
    }

    //obliczenie długości kroku
    double h = (b - a) / static_cast<double>(n);
    //zaalokowanie zmiennej na wynik
    double** result = new double* [n];
    for (int i = 0; i < n; i++) result[i] = new double[2];
    //zmienne pomocnicze do zwiększenia wydajności
    double px = a + h;
    double py;
    double k1 = h * f(a, y0);
    double k2 = h * f(a + (0.5 * h), y0 + (0.5 * k1));
    double k3 = h * f(a + (0.5 * h), y0 + (0.5 * k2));
    double k4 = h * f(px, y0 + k3);
    //obliczenie współrzędnych punktów
    result[0][0] = px;
    result[0][1] = y0 + ((k1 + (2 * k2) + (2 * k3) + k4) / 6.);
    n--;
    for (int i = 1; i < n; i++)
    {
        py = result[i - 1][1];
        double k1 = h * f(px, py);
        double k2 = h * f(px + (0.5 * h), py + (0.5 * k1));
        double k3 = h * f(px + (0.5 * h), py + (0.5 * k2));
        px += h;
        result[i][0] = px;
        double k4 = h * f(px, py + k3);
        result[i][1] = py + ((k1 + (2 * k2) + (2 * k3) + k4) / 6.);
    }
    //przypisanie ,,na sztywno" wartości współrzędnych ostatniego punktu
    result[n][0] = b;
    if (n >= 1)
    {
        py = result[n - 1][1];
        double k1 = h * f(b - h, py);
        double k2 = h * f(b - (0.5 * h), py + (0.5 * k1));
        double k3 = h * f(b - (0.5 * h), py + (0.5 * k2));
        double k4 = h * f(b, py + k3);
        result[n][1] = py + ((k1 + (2 * k2) + (2 * k3) + k4) / 6.);
    }
    //zwrócenie wyniku
    return result;
}

//EIGENVECTOR POWER METHOD
//argumenty: A - wskaźnik na macierz, n - liczba wierszy macierzy A,
//m - liczba liczba kolumn macierzy A, steps - maksymalna liczba kroków
double* eigenvectorPM(double** A, int n, int m, int steps)
{
    //walidacja danych wejściowych
    if (A == nullptr)
    {
        printf("\nWrong arguments! (pointer to matrix)");
        return nullptr;
    }
    if (steps < 1)
    {
        printf("\nWrong arguments! (number of steps)");
        return nullptr;
    }
    if (n < 1)
    {
        printf("\nWrong arguments! (size of matrix)");
        return nullptr;
    }
    if (m < 1)
    {
        printf("\nWrong arguments! (size of matrix)");
        return nullptr;
    }

    //zmienne potrzebne do działań
    //tablica na wektor własny oraz wartość własną - x[m] to miejsce na wartość własną
    double* x = new double[m + 1];
    for (int i = 0; i < m; i += 2) x[i] = 1.;
    for (int i = 1; i < m; i += 2) x[i] = -1.;
    //tablica na wektor pomocniczy
    double* v = new double[n + 1];
    //pomocnicze zmienne typu double
    double value = 0.;
    double sum = 0.;
    //ustawienie parzystej liczby kroków
    if (steps % 2 == 1) steps++;
    //wyliczenie wektora własnego i odpowiadającej mu wartości własnej z zadaną dokładnością
    for (int i = 0; i < steps; i++)
    {
        value = 0.;
        for (int j = 0; j < n; j++)
        {
            //przemnożenie macierzy A przez aktualny wektor x i wpisanie wyniku do wektora v
            sum = 0.;
            for (int k = 0; k < m; k++) sum += A[j][k] * x[k];
            v[j] = sum;
            //wybranie jako wartość własna elementu wektora v z największym modułem
            if (value < abs(sum)) value = abs(sum);
        }
        //przypisanie wartości własnej do odpowiedniej zmiennej
        x[m] = value;
        //warunek stopu
        if (value == 0.)
        {
            //dealokacja pamięci
            delete[] v;
            //zwrócenie wyniku
            return x;
        }
        //przypisanie wektorowi x wartości do kolejnej iteracji
        for (int h = 0; h < m; h++) x[h] = v[h] / value;
    }
    //sprawdzenie znaku wartości własnej
    sum = 0.;
    for (int k = 0; k < m; k++) sum += A[0][k] * x[k];
    double d = x[m] * x[0] * (-1.);
    //dokładność sprawdzenia znaku
    double EPSYLUN = 0.;// 000000001;
    if (abs(sum - d) < EPSYLUN) x[m] = x[m] * (-1.);
    //dealokacja pamięci
    delete[] v;
    //zwrócenie wyniku
    return x;
}