#include <vector>
#include <map>

std::map<double, std::string> element_map = {
    {1.0, "H"},
    {6.0, "C"}
};

// struct bond_type_parameters
// {
    
//     double lambda_ss;
//     double lambda_sp;
//     double lambda_pp_sigma;
//     double lambda_pp_pi;
//     double beta_ss;
//     double beta_sp;
//     double beta_pp_sigma;
//     double beta_pp_pi;
//     double alpha;
//     double gamma;
//     double omega;
//     double rho; 

 
// };
// std::map<std::string, bond_type_parameters> bond_type_parameters_dict = {
//     {"HH", bond_type_parameters{
//         0.280, //lambda_ss
//         0.0, //lambda_sp
//         0.0, //lambda_pp_sigma
//         0.0, //lambda_pp_pi
//         -4.442, //beta_ss
//         0.0, //beta_sp
//         0.0, //beta_pp_sigma
//         0.0, //beta_pp_pi
//         2.823, //alpha
//         12.612, //gamma
//         -0.0791, //omega
//         2.279 //rho
//     }},
//     {"CH", bond_type_parameters{
//         0.275, //lambda_ss
//         0.218, //lambda_sp
//         0.0, //lambda_pp_sigma
//         0.0, //lambda_pp_pi
//         -8.574, //beta_ss
//         -6.813, //beta_sp
//         0.0, //beta_pp_sigma
//         0.0, //beta_pp_pi
//         2.831, //alpha
//         99.370, //gamma
//         -0.0340, //omega
//         2.843 //rho
//     }},
//     {"CC", bond_type_parameters{
//         0.086, //lambda_ss
//         0.180, //lambda_sp
//         0.186, //lambda_pp_sigma
//         0.282, //lambda_pp_pi
//         -5.969, //beta_ss
//         -6.160, //beta_sp
//         -8.420, //beta_pp_sigma
//         -7.403, //beta_pp_pi
//         3.401, //alpha
//         658.659, //gamma
//         0.0312, //omega
//         3.044 //rho
//     }}

// };


std::map<std::string, std::map<std::string, std::pair<double, double>>> bond_type_parameters = {
    {"HH", {
        {"ss", {-4.442, 0.280}}
    }},
    {"CH", {
        {"ss", {-8.574, 0.275}},
        {"sp", {-6.813, 0.218}}
    }},
    {"CC", {
        {"ss", {-5.969, 0.086}},
        {"sp", {-6.160, 0.180}},
        {"pp_sigma", {-8.420, 0.186}},
        {"pp_pi", {-7.403, 0.282}}
    }}
};

std::map<std::string, std::map<std::string, double>> atomic_parameters = {
    {"H", {{"s", -13.605}}},
    {"C", {{"s", -21.559}, {"p", -13.507}}}
};







