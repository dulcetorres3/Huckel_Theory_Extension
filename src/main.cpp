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

struct atom_info
{
    std::string atom;
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
    * @brief Creates a vector of atomic orbitals corresponding to a single atom.
    * 
    * @param atom The name or symbol of the atom (e.g., "H", "C").
    * @param atomic_number The atomic number (Z) of the element.
    * @param coordinates The 3D coordinates of the atom (arma::vec).
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
    /**
    * @brief Creates a vector of bond structures. 
    * 
    * @param aid1 std::vector<double> containing the atom IDs that make up half a bond
    * @param aid2 std::vector<double> containing the atom IDs that make up half a bond
    * @param bond_order std::vector<double> containing the bond orders
    * @param atom_id Atom ID of interest 
    * @return std::vector<bond> A list of bond structures that correspond to the given atom_id; 
    * each bond structures contains the atom-id of the counterpart atom and the bonder order.
    * 
    * 
    */
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
    std::map<double, std::vector<bond>> bond_map; // contains the a vector of all all bond information (atom pair and bond order) corresponding to an atom id
    //std::vector<atomic_orbital> atomic_orbitals;
    double bohr_radius = 0.52917; // units of angstrom
    double k_cal_conversion = 23.06054;
    int total_electrons = 0;
    std::vector<atom_info> atom_vec;
    
    public:
    std::vector<atomic_orbital> atomic_orbitals;
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

            // update total_electrons
            total_electrons+= atomic_number;

            // update atom_vec 
            atom_info my_atom;
            my_atom.atom = atom;
            my_atom.coordinates = coordinates;
            atom_vec.push_back(my_atom);

            // create dictionary entry
            bond_map[atom_id] = get_bond_vector(aid1, aid2, bond_orders, atom_id);
      

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
        /**
        * @brief Determines a double bond as sigma or pi if overlap occurs between orbitals   
        * 
        * @param u A p orbital
        *  @param v A p orbital
        * @return std::string pi bond type (sigma, pi, none-- if no overlap occurs between p orbitals) 
        * 
        * 
        */
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
        /**
        * @brief Determines if a bond occurs between the atoms corresponding to the provided orbitals  
        * 
        * @param u A p orbital
        * @param v A p orbital
        * @return std::tuple<bool, double> bool (true if a bond occurs), double (bond order -- 0.0 if no bond occurs)
        * 
        * 
        */

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
        /**
        * @brief Determines the lambda and beta parameters for the two-center off-diagonal H values  
        * 
        * @param u A p orbital
        * @param v A p orbital
        * @param bond_order the bond order shared between atoms of orbital u and orbital v  
        * @return std::tuple<double, double> containing the lambda and beta values 
        * 
        * 
        */
        std::string u_atom = u.atom;
        std::string v_atom = v.atom;
        std::string u_orbital = u.orbital_type;
        std::string v_orbital = v.orbital_type;
        std::string atom_pair = u_atom + v_atom;
        std::string orbital_pair = u_orbital + v_orbital;
        std::pair<double, double>  bond_parameters;

        if(bond_order==2.0 && orbital_pair=="pp") // determine double bond parameters for pp orbitals 
        {
            std::string bond_type = pi_bond_type(u, v);
            if(bond_type=="none") // "none" means no overlap occured between orbitals (i.e px-py, px-pz, etc.)
            {
                bond_parameters = std::make_pair(0.0, 0.0); 
        
            }
            else{
                bond_parameters = bond_type_parameters[atom_pair][bond_type];
            }
        }
        else if(bond_order==1.0 && orbital_pair=="pp") // determine single bond parameters for pp orbitals 
        {
            if(u.orbital_axis!=v.orbital_axis) // orbitals are along different axis aligns and NO sigma overlap will occurs
            {
                bond_parameters = std::make_pair(0.0, 0.0);
            }
            else // orbitals are along same axis and sigma overlap DOES occur
            {
                orbital_pair = "pp_sigma";
                bond_parameters = bond_type_parameters[atom_pair][orbital_pair]; 
            }
        }
        else // determine single bond parameters (ss, sp)
        {
            bond_parameters = bond_type_parameters[atom_pair][orbital_pair];
        }

        return bond_parameters;
    
    }


    double two_center_off_diagonal(atomic_orbital u, atomic_orbital v)
    {
        /**
        * @brief Determines the two-center off-diagonal H values  between the given orbitals 
        * 
        * @param u A p orbital
        * @param v A p orbital 
        * @return  two-center off-diagonal H value
        * 
        * 
        */

        // determine if orbitals share a bond
        std::tuple<bool, int> bonded = are_atoms_bonded(u,v);
        if(std::get<0>(bonded) == false)
        {
            return 0.0;
        }

        double bond_order = std::get<1>(bonded);
        std::pair<double, double> bond_parameters = determine_bond_parameters(u,v,bond_order);
        double beta = bond_parameters.first;
        double lambda = bond_parameters.second;
        double R_AB = arma::norm(u.coordinates - v.coordinates); // euclidian distance between atoms
        
        return beta * std::pow((R_AB/bohr_radius), 0.5) * std::exp(-lambda * std::pow(R_AB,2) / std::pow(bohr_radius,2));
        
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




// ad this to constructor:
// store coordinates once per atom
atom_coordinates[atom_id] = coordinates;

std::map<double, std::string> atom_element;
atom_element[atom_id] = atom;


double compute_G_AB(std::string A, std::string B, double R_AB) {
    std::string atoms = A + B;
    if (atom_bond_type_parameters.find(atoms) == atom_bond_type_parameters.end()) {
        atoms = B + A;
        if (atom_bond_type_parameters.find(atoms) == atom_bond_type_parameters.end()) {
            return 0.0;
        }
    }

    auto params = atom_bond_type_parameters[atoms];
    double lambda = params[0];
    double gamma = params[1];
    double omega = params[2];
    double r = params[3];

    return gamma * std::exp(-lambda * R_AB)
         + omega * std::exp(-6.0 * std::pow(R_AB - r, 2));
}
double compute_E_rep() {
    double E_rep = 0.0;
    for (auto it1 = atom_coordinates.begin(); it1 != atom_coordinates.end(); ++it1) {
        double id1 = it1->first;
        arma::vec Ri = it1->second;
        std::string Ai = atom_element[id1];

        auto it2 = it1;
        ++it2;
        for (; it2 != atom_coordinates.end(); ++it2) {
            double id2 = it2->first;
            arma::vec R_j = it2->second;
            std::string A_j = atom_element[id2];

            double R_AB = arma::norm(Ri - Rj);
            E_rep += compute_G_AB(Ai, Aj, R_AB);
        }
    }
    return E_rep;
}



double compute_total_energy() {
    arma::mat Hmat = H();
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, Hmat);

    int num_electrons = atomic_orbitals.size(); 
    int num_occ_orbitals = num_electrons / 2;

    double E_electronic = arma::sum(eigval.subvec(0, num_occ_orbitals - 1));
    double E_rep = compute_E_rep();

    return E_electronic + E_rep;
}





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




