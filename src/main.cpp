#include <nlohmann/json.hpp> 
#include <iostream>
#include <vector>
#include <fstream>
#include <armadillo>
#include <filesystem>
#include "basis_sets.h"
#include <tuple>

using json = nlohmann::json;

struct atomic_orbital
{
    double atom_id;
    std::string atom;
    std::string orbital_type;
    std::string orbital_axis;
    double atomic_number;
    arma::vec coordinates;

};

struct bond
{
    double atom_id; // the atom_id of the counterpart
    double bond_order;
};

std::vector<atomic_orbital> create_atomic_orbitals(int atom_id, std::string atom, double atomic_number, arma::vec coordinates, bool hydrogen=false)
{
    /**
    * @brief Creates a vector of atomic orbitals based on input parameters.
    * 
    * @param atom_ The name or symbol of the atom (e.g., "H", "C").
    * @param atomic_number_ The atomic number (Z) of the element.
    * @param coordinates_ The 3D coordinates of the atom (arma::vec).
    * @return std::vector<atomic_orbital> A list of atomic orbitals.
    * 
    * 
    */

    std::vector<atomic_orbital> atomic_orbitals;

    atomic_orbital s_orbital;
    s_orbital.atom_id = atom_id;
    s_orbital.atom = atom;
    s_orbital.orbital_type = "s";
    s_orbital.atomic_number = atomic_number;
    s_orbital.coordinates = coordinates;
    s_orbital.orbital_axis = "none"; // s orbitals are not along any axis
    atomic_orbitals.push_back(s_orbital);

    if(!hydrogen)
    {
        for(int i=0; i < 3; i++)
        {
            atomic_orbital p_orbital;
            p_orbital.atom_id = atom_id;
            p_orbital.atom = atom;
            p_orbital.orbital_type = "p";
            p_orbital.atomic_number = atomic_number;
            p_orbital.coordinates = coordinates;
            p_orbital.orbital_axis = orbital_axis_map[i + 1]; // x axis= 1, y axis = 2, z axis = 3
            atomic_orbitals.push_back(p_orbital);

        }
    }
   
    return atomic_orbitals;
}

std::vector<bond> get_bond_vector(std::vector<double> aid1, std::vector<double> aid2, std::vector<double> bond_order, double atom_id)
{
    std::vector<bond> my_bonds;

    for(int i=0; i < aid1.size(); i++)
    {
        double atom_id1 = aid1[i];
        double atom_id2 = aid2[i];
        if(atom_id1==atom_id)
        {
            bond bond_counterpart;
            bond_counterpart.atom_id = aid2[i];
            bond_counterpart.bond_order = bond_order[i];
            my_bonds.push_back(bond_counterpart);
        }
        else if(atom_id2==atom_id)
        {
            bond bond_counterpart;
            bond_counterpart.atom_id = aid1[i];
            bond_counterpart.bond_order = bond_order[i];
            my_bonds.push_back(bond_counterpart);
        }
    }
    return my_bonds;
}

class Extended_Huckel
{
    private:
    int N;
    std::map<double, std::vector<bond>> bond_map; 
    std::vector<atomic_orbital> atomic_orbitals;
    double sigma_threshold = 0.9;
    double pi_threshold = 0.1;
    double bohr_radius = 0.52917;
    
    public:
    //std::vector<atomic_orbital> atomic_orbitals;
    Extended_Huckel(std::vector<std::vector<double>> compound_data)
    {
        //coordinate information
        std::vector<double> atomic_numbers = compound_data[0];
        std::vector<double> x_coords = compound_data[1];
        std::vector<double> y_coords = compound_data[2];
        std::vector<double> z_coords = compound_data[3];
        //bond information
        std::vector<double> aid1 = compound_data[4];
        std::vector<double> aid2 = compound_data[5];
        std::vector<double> bond_orders = compound_data[6];

        int total_atoms = atomic_numbers.size();

       
        for(int i=0; i < total_atoms; i++)
        {
            double atomic_number = atomic_numbers[i];
            std::string atom = element_map[atomic_number];
            double atom_id = i + 1.0;
            arma::vec coordinates = {x_coords[i], y_coords[i], z_coords[i]};

            // create dictionary entry
            bond_map[atom_id] = get_bond_vector(aid1, aid2, bond_orders, atom_id);
            std::cout << "atom id: " << atom_id << std::endl;
            std::vector<bond>  my_bonds = bond_map[atom_id];
            std::cout << "MY BONDS: " << std::endl;
            for(int j=0; j< my_bonds.size(); j++)
            {
                std::cout << my_bonds[j].bond_order << std::endl;
            }

            // create the atomic orbitals for atom
            if(atom=="H")
            {
                std::vector<atomic_orbital> my_orbitals = create_atomic_orbitals(atom_id, atom, atomic_number, coordinates, true);
                for (const auto& atomic_orbital : my_orbitals) {
                    atomic_orbitals.push_back(atomic_orbital);
                }
            }
            else
            {
                std::vector<atomic_orbital> my_orbitals = create_atomic_orbitals(atom_id, atom, atomic_number, coordinates);
                for (const auto& atomic_orbital : my_orbitals) {
                    atomic_orbitals.push_back(atomic_orbital);
                }
            }
        }

        N = atomic_orbitals.size();


    }

    
    std::string pi_bond_type(atomic_orbital u, atomic_orbital v)
    {
        std::string u_orbital_axis = u.orbital_axis;
        std::string v_orbital_axis = v.orbital_axis;
        std::string bond_type = "none"; 

        if(u_orbital_axis!=v_orbital_axis) // mixed orbital interactions will have no overlap (i.e. px-py, pz-py, etc)
        {
            return bond_type;
        }
        else if(u_orbital_axis==v_orbital_axis && u_orbital_axis=="z") // matching orbitals along the z axis are labeled as sigma 
        {
            bond_type = "pp_sigma";
        }
        else if(u_orbital_axis==v_orbital_axis && (u_orbital_axis=="x" || u_orbital_axis=="y")) // matching orbital along either x or y axis are labeled as pi
        {
            bond_type = "pp_pi";
        }

        return bond_type;
    }

    std::tuple<bool, double> are_atoms_bonded(atomic_orbital u, atomic_orbital v)
    {
        double A_atom_id = u.atom_id;
        double B_atom_id = v.atom_id; 
        bool bonded = false;
        double bond_order = 0.0;
        std::vector<bond> A_bonds = bond_map[A_atom_id];

        for(int i=0; i < A_bonds.size(); i++)
        {
            double counterpart_id = A_bonds[i].atom_id;
            if(counterpart_id == B_atom_id) 
            {
                bonded = true;
                bond_order = A_bonds[i].bond_order;
            }
        }
        std::tuple<bool, double> bond_elements(bonded, bond_order);
        return bond_elements;
    }

    std::pair<double, double> determine_bond_parameters(atomic_orbital u, atomic_orbital v, double bond_order)
    {
        std::string u_atom = u.atom;
        std::string v_atom = v.atom;
        std::string u_orbital = u.orbital_type;
        std::string v_orbital = u.orbital_type;
        std::string atom_pair = u_atom + v_atom;
        std::pair<double, double>  bond_parameters;
    

        if(bond_order==2.0) // determine double bond parameters
        {
            std::string bond_type = pi_bond_type(u, v);
            if(bond_type=="none") // "none" means no overlap occured between orbitals (i.e px-py, px-pz, etc.)
            {
                bond_parameters = std::make_pair(0.0, 0.0); // 0 values will ensure the off diagonal is 0
        
            }
            else{
                bond_parameters = bond_type_parameters[atom_pair][bond_type];
            }
        }
        else // determine single bond parameters
        {
            bond_parameters = bond_type_parameters[atom_pair][u_orbital+v_orbital];
        }

        return bond_parameters;
    
    }


    double two_center_off_diagonal(atomic_orbital u, atomic_orbital v)
    {
        // determine if orbitals share a bond
        std::tuple<bool, int> bonded = are_atoms_bonded(u,v);
        if(std::get<0>(bonded) == false)
        {
            return 0.0;
        }

        double bond_order = std::get<1>(bonded);
        std::pair<double, double> bond_parameters = determine_bond_parameters(u,v,bond_order);
        double lambda = bond_parameters.first;
        double beta = bond_parameters.second;
        double R_AB = arma::norm(u.coordinates - v.coordinates); // euclidian distance between atoms
        

        return beta * std::pow((R_AB/bohr_radius), 0.5) * std::exp((-lambda * std::pow(R_AB,2)) / std::pow(bohr_radius,2));
        
    }


    arma::mat H()
    {
        arma::mat H_matrix(N,N);

        for(int u=0; u < N; u++)
        {
            for(int v=0; v < N; v++)
            {
                atomic_orbital u_orbital = atomic_orbitals[u];
                atomic_orbital v_orbital = atomic_orbitals[v];
                double u_atom_id = u_orbital.atom_id;
                double v_atom_id = v_orbital.atom_id;
                std::string u_atom = u_orbital.atom;

                // diagonal terms 
                if(u==v)
                {
                    std::string orbital_type = u_orbital.orbital_type;
                    H_matrix(u,v) = atomic_parameters[u_atom][orbital_type];
                }

                // one-center off diagonal terms
                else if(u_atom_id==v_atom_id)
                {
                    H_matrix(u,v) = 0.0;
                }
                // two-center off diagon terms
                else{
                    H_matrix(u,v) = two_center_off_diagonal(u_orbital, v_orbital);
                }
            }
        }
        return H_matrix;
    }


};

std::vector<std::vector<double>> read_json(std::string filepath)
{
    std::ifstream input(filepath);
    if (!input.is_open()){
        std::cerr << "Error: Could not open JSON file!" << std::endl;
        return {};
    }

    json data;
    
    try{
        input >> data;
    } catch (json::parse_error& e){
        std::cerr << "Parse error: " << e.what() << std::endl;
    }
    
    // Extract elements, coordinates, and bond information 
    std::vector<double> elements, x_coords, y_coords, z_coords, aid1, aid2, bond_orders;
    
    // Get the first compound (assuming single compound in array)
    auto& compound = data["PC_Compounds"][0];
    
    // Extract atomic elements
    for (auto& element : compound["atoms"]["element"]) {
        elements.push_back(element);
    }
    
    // Extract coordinates (assuming first conformer)
    auto& conformer = compound["coords"][0]["conformers"][0];
    for (auto& x : conformer["x"]) {
        x_coords.push_back(x);
    }
    for (auto& y : conformer["y"]) {
        y_coords.push_back(y);
    }
    for (auto& z : conformer["z"]) {
        z_coords.push_back(z);
    }
    // Extract bond information
    for (auto& a1 : compound["bonds"]["aid1"]) {
        aid1.push_back(a1);
    }
    for (auto& a2 : compound["bonds"]["aid2"]) {
         aid2.push_back(a2);
    }
    for (auto& order : compound["bonds"]["order"]) {
        bond_orders.push_back(order);
    }

    std::vector<std::vector<double>>  compound_data = {
        elements, // index 0
        x_coords, // index 1
        y_coords, // index 2
        z_coords, // index 3
        aid1, // index 4
        aid2, //index 5
        bond_orders //index 6
    };
    return compound_data;

}

int main(int argc, char** argv)
{

    std::string filepath = argv[1];
    // Read JSON file
    std::vector<std::vector<double>>  compound_data = read_json(filepath);
    Extended_Huckel compound(compound_data);
    std::cout << compound.H() << std::endl;
    
    /*
    int N = compound.atomic_orbitals.size();
    for(int i=0; i < N; i++)
    {
        atomic_orbital orbital = compound.atomic_orbitals[i];
        std::cout << orbital.atom_id << std::endl;
    }
    */
   
    return 0;
}




