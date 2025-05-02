#include <vector>
#include <map>
#include <armadillo>

std::map<double, std::string> element_map = {
    {1.0, "H"},
    {6.0, "C"}
};

std::map<std::string, std::map<std::string, std::pair<double, double>>> bond_type_parameters = {
    {"HH", {
        {"ss", {-4.442, 0.280}} // beta_ss , lambda_ss
    }},
    {"CH", {
        {"ss", {-8.574, 0.275}}, // beta_ss , lambda_ss
        {"sp", {-6.813, 0.218}} // beta_sp , lambda_sp
    }},
    {"CC", {
        {"ss", {-5.969, 0.086}}, // beta_ss , lambda_ss
        {"sp", {-6.160, 0.180}}, // beta_sp , lambda_sp
        {"pp_sigma", {-8.420, 0.186}}, // beta_pp_sigma , lambda_pp_sigma
        {"pp_pi", {-7.403, 0.282}} // beta_pp_sigma , lambda_pp_sigma
    }}
};

std::map<std::string, std::vector<double>> atom_bond_type_parameters = {
    {"HH", {2.823, 12.612, -0.0791, 2.279}}, // alpha, gamma, omerga, r
    {"CH", {2.831, 99.370, -0.0340, 2.843}},
    {"CC", {3.401, 658.659, 0.0312, 3.044}}
};

std::map<std::string, std::map<std::string, double>> atomic_parameters = {
    {"H", {{"s", -13.605}}},
    {"C", {{"s", -21.559}, {"p", -13.507}}},
};

std::map<int, std::string> orbital_axis_map = {
    {1, "x"},
    {2, "y"},
    {3, "z"}
};


std::map<std::string, arma::vec> orbital_axis_coordinates = {
    {"x", arma::vec{1.0, 0.0, 0.0}},  // x-axis unit vector
    {"y", arma::vec{0.0, 1.0, 0.0}},  // y-axis unit vector
    {"z", arma::vec{0.0, 0.0, 1.0}}   // z-axis unit vector
};

