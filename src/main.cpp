#include <iostream>
#include <vector>
#include <fstream>
#include <nlohmann/json.hpp> 
#include <armadillo>
#include <filesystem>
#include "basis_sets.h"

using json = nlohmann::json;

struct atomic_orbital
{
    double atom_id;
    std::string atom;
    std::string orbital_type;
    double atomic_number;
    arma::vec coordinates;
    arma::vec p_direction; //  for p orbitals


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
        std::vector<double> bond_orders = compound_data[5];

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
                std::cout << my_bonds[j].atom_id << std::endl;
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










/**
 * @brief Constructs the Extended Hückel Hamiltonian matrix H (N x N) using the MTB/2 method.
 * 
 * The diagonal elements are assigned from the empirical atomic parameters table (onsite energies).
 * The off-diagonal elements between bonded orbitals are calculated using Voityuk's MTB/2 formula:
 * 
 *     H_{μν} = β * sqrt(R_AB / a₀) * exp(-λ * (R_AB / a₀)^2)
 * 
 * where β and λ are empirical parameters specific to the atom pair and orbital interaction type.
 * 
 * @return arma::mat - the NxN Hamiltonian matrix for the molecular system.
 */

 
arma::mat H_matrix()
{
    arma::mat H(N, N, arma::fill::zeros);


    for (int mu = 0; mu < N; mu++)
    {
        atomic_orbital mu_orbital = atomic_orbitals[mu];
        double mu_atom_id = mu_orbital.atom_id;

        for (int nu = 0; nu < N; nu++)
        {
            atomic_orbital nu_orbital = atomic_orbitals[nu];
            double nu_atom_id = nu_orbital.atom_id;


            // get the diagonal term form the atomic_parameters dictionary for H and C
            if (mu == nu)
{
    std::string atom = mu_orbital.atom;
    std::string orb = mu_orbital.orbital_type;

    if (atomic_parameters.count(atom) && atomic_parameters[atom].count(orb))
    {
        H(mu, nu) = atomic_parameters[atom][orb];
    }
    else
    {
        H(mu, nu) = 0.0; //if the onsite energy atomic paramter is missing
    }
    continue;
}

// for off diag, if s-p or p-p are on the same atom
if (mu_atom_id == nu_atom_id)
{
    H(mu, nu) = 0.0;
    continue;
}

//if atoms are not bonded
if (!bonded_flag)
{
    //  non-bonded elements to zero 
    H(mu, nu) = 0.0;
    H(nu, mu) = 0.0;
    continue;
}

//find the bond type and orbital interaction as ss, sp, pp_sigma, pp_pi
            std::string atom1 = mu_orbital.atom;
            std::string atom2 = nu_orbital.atom;
            std::string bond_type = (atom1 < atom2) ? atom1 + atom2 : atom2 + atom1;

            std::string orb1 = mu_orbital.orbital_type;
            std::string orb2 = nu_orbital.orbital_type;
            std::string orbital_interaction;

            if (orb1 == "s" && orb2 == "s")
            {
                orbital_interaction = "ss";
            }
            else if ((orb1 == "s" && orb2 == "p") || (orb1 == "p" && orb2 == "s"))
            {
                orbital_interaction = "sp";
            }
            else if (orb1 == "p" && orb2 == "p")
            {

//NEED TO differentiate between pp sigma and pp pi...How??
                {
                    orbital_interaction = "pp_sigma";
                }
                else
                {
                    orbital_interaction = "pp_pi";
                }
            }
            else
            {
                continue; 
            }

            //get the beta and lambda parameters using the bond_type_parameters dictionary
            auto bond_params = bond_type_parameters[bond_type][orbital_interaction];
            double beta = bond_params.first;
            double lambda = bond_params.second;

            //solve using equation 5
            arma::vec R_mu = mu_orbital.coordinates;
            arma::vec R_nu = nu_orbital.coordinates;
            double R_AB = arma::norm(R_mu - R_nu);

            const double a0 = 0.52917;

            double H_mu_nu = beta * std::sqrt(R_AB / a0) * std::exp(-lambda * std::pow(R_AB / a0, 2));

            H(mu, nu) = H_mu_nu;
            H(nu, mu) = H_mu_nu;
        }
    }

    return H;
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




