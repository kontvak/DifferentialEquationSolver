#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <vector>
#include <cmath>
#include "json.hpp"
using json = nlohmann::json;


template <int N, typename T>
class State {
private:
    std::array<T, N> x;
public:

    State() : x() {}
    State(std::array<T, N> X) : x(std::move(X)) {}

    ~State() = default;
    State(const State& other) = default;
    State(State&& other) : x(std::move(other.x)) {}

    State& operator=(const State& other) = default;
    State& operator=(State&& other) {
        if (this != &other) {
            x = std::move(other.x);
        }
        return *this;
    }

    State<N, T> operator+(const State<N, T>& other) const {
        std::array<T, N> result;
        for (int i = 0; i < N; ++i) {
            result[i] = x[i] + other.x[i];
        }
        return State<N, T>(result);
    }

    State<N, T> operator*(double scalar) const {
        std::array<T, N> result;
        for (int i = 0; i < N; ++i) {
            result[i] = static_cast<T>(x[i] * scalar);
        }
        return State<N, T>(result);
    }

    //T operator[](int i) { return x[i]; } not working
    const T operator[](int i) const { return x[i]; }
    T& operator[](int i) { return x[i]; }

    void print() {
        for (int i = 0; i < N; ++i) {
            std::cout << x[i] << ' ';
        }
        std::cout << std::endl;
    }
    void print_array(double* x0_ptr, double* x1_ptr) //TODO write for n-dim or write out of class
    {
        *x0_ptr = x[0];
        *x1_ptr = x[1];
    }
};

template <typename state>
class Function {
private:
    std::vector<double> param;
public:
    Function(std::vector<double> param) : param(param) {}
    Function() : param(std::vector<double>()) {}

    state dfdt(double t, const state& a) {
        state derivative_a;
        derivative_a[1] = -param[0] * param[0] * a[0] - param[1]* a[1]+ std::sin(0.5*param[0]*t);
        derivative_a[0] = a[1];
        return derivative_a;
    }
};



template <typename state>
class Method1 {
public:
    Method1() {}
    void make_step(state& s_ptr, Function<state> f, double t, double dt) {
        s_ptr = s_ptr + f.dfdt(t, s_ptr) * dt;
    }
};
template <typename state>
class Method2 {
public:
    Method2() {}
    void make_step(state& s_ptr, Function<state> f, double& t, double dt) {
        state tmp = s_ptr;
        state delta = f.dfdt(t, tmp)*dt;

        s_ptr = s_ptr + (f.dfdt(t, tmp + delta) * dt + delta) * 0.5;
        t += dt;
    }
};

/*
template <typename state>
void make_step(state& s_ptr, Function<state> f, double t, double dt) {
    s_ptr = s_ptr + f.dfdt(t, s_ptr) * dt;
}
*/



int main(int argc, char* argv[]) {

    std::ifstream f(argv[1]);
    json json_data = json::parse(f);
    f.close();

    int N = json_data["N"];
    double STEP = json_data["STEP"];
    double W = json_data["W"]; //TODO more parametrs from config
    std::string method = json_data["method"];
    std::string out_file = json_data["out_file"];

    std::array<double, 2> b{ 1., 0. };
    State<2, double> st(b);

    /*TODO write class Abstract method
    if (method == std::string("1")) {
        Method1<State<2, double>> meth;
    }
    else if (method == std::string("2")) {
        Method2<State<2, double>> meth;
    }
    */
    Method2<State<2, double>> meth;
    Function<State<2, double>> fu({ W , 0.1});
    double t = 0;


    double* t_arr{ new double[N] {} };
    double* x0_arr{ new double[N] {} };
    double* x1_arr{ new double[N] {} };

    for (auto i = 0; i < N; ++i) {

        st.print_array( x0_arr + i, x1_arr + i);
        *(t_arr + i) = t;
        meth.make_step(st, fu, t, STEP);
    }

    std::ofstream tabel1;
    tabel1.open(out_file);
    tabel1 << "t" << ',' << "x0" << ',' << "x1" << '\n';
    for (auto i = 0; i < N; ++i) {
        tabel1 << t_arr[i] << ',' << x0_arr[i] << ',' << x1_arr[i] << '\n';
    }
    tabel1.close();
    delete[] t_arr;
    delete[] x0_arr;
    delete[] x1_arr;
}
