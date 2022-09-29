/* Just postprocessing of several RDFs in text-files. */

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <sstream>
#include <cmath>
#include <iterator>
#include <memory>
#include <algorithm>
#include <string>


typedef std::tuple<double, double> data_tuple;


std::vector <double> mesh (double left_border, const double & right_border, const double & step);

std::vector<double> cubic_spline (std::vector<data_tuple> & data, std::vector<double> & xx);

std::string exec (const std::string& str);

template<typename T>
T fromString (const std::string& s);

std::vector<data_tuple> coordinates_read (const std::string & name);

std::vector<data_tuple> collection (const std::string & files_name, int & data_files_count);

void averaging (std::vector<data_tuple> & data, int & data_count);

void data_file_creation (const std::string & name, std::vector<double> & x, std::vector<double> & y);

void plot (const std::string & name, const double & left, const double & right,
           const std::string & title, const std::string & xlabel, const std::string & ylabel);


int main () {
    int data_files_count = fromString<int>(exec("find ./gOH -type f | wc -l"));
    std::vector<data_tuple> collected_data = std::move(collection("gOH/gOH.", data_files_count));
    averaging(collected_data, data_files_count);

    double left_border = std::get<0>(collected_data[0]);
    double right_border = std::get<0>(collected_data[collected_data.size()-1]);
    double mesh_step = (right_border-std::get<0>(collected_data[collected_data.size()-2]))/10.00;

    std::vector<double> xx = std::move(mesh(left_border, right_border, mesh_step));
    std::vector<double> yy = std::move(cubic_spline(collected_data, xx));
    data_file_creation("splined.txt", xx, yy);
    plot ("splined", left_border, right_border, "RDF", "r, Angstroms", "a.u.");
    return 0;
}


//The function returns the terminal ans. Input - string for term.
std::string exec (const std::string& str) {
    const char* cmd = str.c_str();
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    result = result.substr(0, result.length()-1);
    return result;
}


// Returns number from string.
template<typename T>
T fromString (const std::string& s) {
    std::istringstream iss(s);
    T res;
    iss >> res;
    return res;
}


namespace std {
    istream& operator >> (istream& in, data_tuple & coordinates) {
        double first, second;
        in >> first >> second;
        coordinates = {first, second};
        return in;
    }

    ostream& operator << (ostream& out, const data_tuple & coordinates) {
        auto [first, second] = coordinates;
        out << first << ' ' << second << ' ';
        return out;
    }
}

// Read data from columns in text-file.
std::vector<data_tuple> coordinates_read (const std::string & name) {
    std::ifstream fin(name);
    if (!fin.is_open()) throw std::runtime_error("Error opening file.");
    std::vector<data_tuple> tuples_vector;
    copy(std::istream_iterator<data_tuple> {fin},
         std::istream_iterator<data_tuple> {},
         back_inserter(tuples_vector));
    //copy(tuples_vector.begin(), tuples_vector.end(), std::ostream_iterator<data>(std::cout, "\n"));
    return tuples_vector;
}


// std::to_string not safe enough. It will be used everywhere instead of std::to_string.
template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


// Collects data in tuples for averaging.
template<size_t Is = 0, typename... Tp>
void sum_coordinates (std::tuple<Tp...>& coordinate, std::tuple<Tp...>& new_data) {
    std::get<Is>(coordinate) += std::get<Is>(new_data);
    if constexpr(Is + 1 != sizeof...(Tp))
        sum_coordinates <Is + 1>(coordinate, new_data);
}


// Collects data from all files in one std::vector<data_tuple>.
std::vector<data_tuple> collection (const std::string & files_name, int & data_files_count) {
    std::vector<data_tuple> data = coordinates_read(files_name + toString(0));
    if (data_files_count == 1) return data;
    for (int i = 1; i < data_files_count; ++i)
        for (int j = 0; j < data.size(); ++j) {
            std::vector<data_tuple> new_data = std::move(coordinates_read(files_name + toString(i)));
            sum_coordinates(data[j], new_data[j]);
        }
    return data;
}


// Averaging tuples content.
template<size_t Is = 0, typename... Tp>
void averaged_coordinates (std::tuple<Tp...>& coordinate, const double & data_count) {
    std::get<Is>(coordinate) /= data_count;
    if constexpr(Is + 1 != sizeof...(Tp))
        averaged_coordinates <Is + 1>(coordinate, data_count);
}

// Averaging data in std::vector<data_tuple>.
void averaging (std::vector<data_tuple> & data, int & data_count) {
    for (auto & i : data)
        averaged_coordinates(i, data_count);
}


// Creates std::vector<double> of points between left_border and right_border with given step.
std::vector <double> mesh (double left_border, const double & right_border, const double & step) {
    std::vector <double> xx ((right_border-left_border) / step);
    std::generate(xx.begin(), xx.end(), [&] {left_border += step; return left_border;});
    return xx;
}


// Creates two std::vector<double> from std::vector<data_tuple>.
void data_from_tuple_vector (std::vector<data_tuple> & data, std::vector<double> & x, std::vector<double> & y) {
    for (auto & i : data) {
        x.emplace_back(std::get<0>(i));
        y.emplace_back(std::get<1>(i));
    }
}


// Two functions below need for creation of cubic spline (like matlab spline) for representation.
std::vector<double> spline_moments (std::vector<double> & f, double & h) {
    std::vector<double> m;
    int n = f.size()-1;
    m.emplace_back((3*f[0] - 4*f[1] + f[2])/2.0/h);
    for (int i = 1; i < n; ++i)
        m.emplace_back((f[i+1] - f[i-1])/2.0/h);
    m.emplace_back((-3*f[n] + 4*f[n-1] - f[n-2])/2.0/h);
    return m;
}

std::vector<double> cubic_spline (std::vector<data_tuple> & data, std::vector<double> & xx) {
    std::vector<double> x, f, yy;
    data_from_tuple_vector(data, x, f);
    int i = 0;
    double h = x[1] - x[0];
    std::vector<double> m = std::move(spline_moments(f, h));
    for (int j = 0; j < x.size()-1; ++j)
        while (xx[i] >= x[j] && xx[i] <= x[j+1]) {
            double buf1 = x[j+1] - xx[i];
            double buf2 = xx[i] - x[j];
            yy.emplace_back(pow(buf1, 2)*(2*buf2+h)*f[j]/pow(h, 3)+
                            pow(buf2, 2)*(2*buf1+h)*f[j+1]/pow(h, 3)+
                            pow(buf1, 2)*buf2*m[j]/pow(h, 2)-
                            pow(buf2, 2)*buf1*m[j+1]/pow(h, 2));
            ++i;
        }
    return yy;
}


// Creates data-file from two std::vector<double> of coordinates with given name.
void data_file_creation (const std::string & name, std::vector<double> & x, std::vector<double> & y) {
    std::ofstream fout;
    fout.open(name, std::ios::trunc);
    for (int i = 0; i < x.size(); ++i)
        fout << toString(x[i]) << '\t' << toString(y[i]) << '\n';
    fout.close();
}


// Creates plot from text-file of coordinates via GNUPlot.
void plot (const std::string & name, const double & left, const double & right,
           const std::string & title, const std::string & xlabel, const std::string & ylabel) {
    std::string range = "[" + toString(left) + ":" + toString(right) + "]";
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 700, 700",
                                      "set output \'" + name + ".jpg\'",
                                      "set title \'" + title + "\'",
                                      "set grid xtics ytics",
                                      "set xrange " + range,
                                      "set xlabel \'" + xlabel + "\'",
                                      "set ylabel \'" + ylabel + "\'",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "plot \'" + name + ".txt\' using 1:2 w lines",
                                      "set terminal pop",
                                      "set output",
                                      "replot", "q"};
    for (const auto & it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}
