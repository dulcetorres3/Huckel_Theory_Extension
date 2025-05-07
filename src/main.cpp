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
    double k_cal_conversion = 23.06054; // converts eV to kcal/mol
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

        // assign the orbitals for each atom 
        for(int i=0; i < total_atoms; i++)
        {
            double atomic_number = atomic_numbers[i];
            std::string atom = element_map[atomic_number];
            double atom_id = i + 1.0;
            arma::vec coordinates = {x_coords[i], y_coords[i], z_coords[i]};

            // update total_electrons
            if(atomic_number==6.0)
            {
                atomic_number=4.0; // accounts for valence C electrons
            }
            total_electrons+= atomic_number;

            // update atom_vec 
            atom_info my_atom;
            my_atom.atom = atom;
            my_atom.coordinates = coordinates;
            atom_vec.push_back(my_atom);

            // create dictionary entry
            bond_map[atom_id] = get_bond_vector(aid1, aid2, bond_orders, atom_id); // entry contains all bonds the atom is involved in 
      

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

    
    std::string pp_bond_type(atomic_orbital u, atomic_orbital v)
    {
        /**
        * @brief Determines a double bond as sigma or pi   
        * 
        * @param u A p orbital
        * @param v A p orbital
        * @return std::string pi bond type (sigma, pi, none-- if no overlap occurs between p orbitals) 
        * 
        * 
        */
        std::string u_orbital_axis = u.orbital_axis;
        std::string v_orbital_axis = v.orbital_axis;
        std::string orbital_pair = u.orbital_type + v.orbital_type;
        //std::cout << "orbital pair: " << orbital_pair << std::endl;
        std::string bond_type;

        if(orbital_pair!="pp") // acounts for sp/ps interactions
        {
            bond_type = orbital_pair;
        }
        else if(u_orbital_axis!=v_orbital_axis) // accounts for mixed pp orbital interactions
        {
            bond_type = "pp_pi";
        }
        else // acounts for identicaly pp orbital interactions
        {
            arma::vec coordinate_difference = u.coordinates - v.coordinates;
            double magnitude = arma::norm(coordinate_difference);
            arma::vec unit_vec = coordinate_difference/magnitude;
            int my_index = orbital_index[u_orbital_axis];
            double unit_vec_component = unit_vec(my_index);
            if(unit_vec_component > 0.9)
            {
                bond_type = "pp_sigma";
            }
            else{
                bond_type = "pp_pi";
            }
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

        if(orbital_pair=="pp") // determine bond parameters for pp orbitals 
        {
            std::string bond_type = pp_bond_type(u, v);
            bond_parameters = bond_type_parameters[atom_pair][bond_type];
            
        }
        else // determine bond parameters for ss/sp orbitals
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
        /*
        if(std::get<0>(bonded) == false)
        {
            return 0.0;
        }
        */
        double bond_order = std::get<1>(bonded);
        std::pair<double, double> bond_parameters = determine_bond_parameters(u,v,bond_order);
        double beta = bond_parameters.first;
        double lambda = bond_parameters.second;
        //std::cout << "beta: " << beta << std::endl;
        //std::cout << "lambda: " << lambda << std::endl;
        double R_AB = arma::norm(u.coordinates - v.coordinates); // euclidian distance between atoms in angstrom 
        //std::cout << "R_AB: " << R_AB  << std::endl;
        //std::cout << "first term: " << beta * std::pow((R_AB/bohr_radius), 0.5) << std::endl;
        //std::cout << "second term: " << std::exp(-lambda * std::pow(R_AB,2) / std::pow(bohr_radius,2)) << std::endl;
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


    double occupation_level(int level)
    {
        int double_occupied_threshold = total_electrons / 2;

        if(level < double_occupied_threshold)
        {
            return 2.0;
        }
        else if(level == double_occupied_threshold && (total_electrons % 2 != 0 ))
        {
            return 1.0;
        }
        else{
            return 0.0;
        }
    }

    double energy_summation()
    {
        arma::mat one_electron_hamiltonian = H();

        arma::mat C;
        arma::vec energies;
        double total_energy = 0.0;

        // solve for energies
        arma::eig_sym(energies, C, one_electron_hamiltonian);
        
        for(int i=0; i < N; i++)
        {
            double energy = energies(i);
            double n = occupation_level(i);
            total_energy+= n*energy;
        }
        return total_energy;

    }

    double G_AB(std::string A, std::string B, double R_AB)
    {
        std::string atom_pair = A + B;
        std::vector<double> parameters = atom_bond_type_parameters[atom_pair]; // alpha, gamma, omerga, r
        double alpha = parameters[0];
        double gamma = parameters[1];
        double omega = parameters[2];
        double r = parameters[3];

        return (gamma * std::exp(-alpha*R_AB)) + (omega * std::exp(-6 * std::pow(R_AB-r,2)));
    }

    double repulsion_energy()
    {
        double summation = 0.0;
        for(int i=0; i < atom_vec.size(); i++)
        {
            for(int j= i+1; j < atom_vec.size(); j++)
            {
                atom_info A_info = atom_vec[i];
                atom_info B_info = atom_vec[j];

                std::string A = A_info.atom;
                std::string B = B_info.atom;
                double R_AB = arma::norm(A_info.coordinates - B_info.coordinates);
            
                summation += G_AB(A, B, R_AB);
            }
        }
        return summation;
    }

    double total_energy()
    {
        double orbital_energy = energy_summation();
        double repulsion = repulsion_energy();
        //std::cout << "orbital energy: " << orbital_energy << std::endl;
        //std::cout << "repulsion: " << repulsion << std::endl;
        //std::cout << "total energy: " << orbital_energy + repulsion << std::endl;
        return (orbital_energy + repulsion) * k_cal_conversion;
    }

    double E_isol_summ()
    {
        double sum = 0.0;

        for(int i =0; i < atom_vec.size(); i++)
        {
            atom_info A_info = atom_vec[i];
            std::string A = A_info.atom;
            std::vector<double> electrons = valence_electrons[A];
            double n_s = electrons[0];
            double n_p = electrons[1];

            double U_s = atomic_parameters[A]["s"];
            double U_p = 0.0;
            if(A == "C")
            {
                U_p = atomic_parameters[A]["p"];
            }

            sum += (n_s * U_s) + (n_p * U_p);
            
        }

        return sum;
    }

    double experimental_heat_summ()
    {
        double summ=0.0;
        for(int i =0; i < atom_vec.size(); i++)
        {
            atom_info A_info = atom_vec[i];
            std::string A = A_info.atom;
            summ+= heat_formation_map[A];
        }
        return summ;
    }

    double enthalpy()
    {
        double energy = k_cal_conversion * total_energy();
        double E_isol = k_cal_conversion * E_isol_summ();
        double heat_formation = experimental_heat_summ();
        //std::cout << "energy (ev): " << total_energy() << std::endl;
        //std::cout << "energy (kcal/mol): " << energy << std::endl;
        //std::cout << "isolation energy (kcal/mol): " << E_isol << std::endl;
        //std::cout << "heat formation (kcal/mol): " << heat_formation << std::endl;
        //std::cout << "enthalpy(kcal/mol): " << energy - E_isol + heat_formation << std::endl;
        return energy - E_isol + heat_formation;
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
    compound.enthalpy();
    // hydrogen data (pubchem did not have downloadable 3d coordinates for hydrogen)
    
    std::vector<double> atomic_numbers = {1, 1};
    std::vector<double> x_coords = {0.0, 1.3984};
    std::vector<double> y_coords = {0.0, 0.0};
    std::vector<double> z_coords = {0.0, 0.0};
    std::vector<double> aid1 = {1.0};
    std::vector<double> aid2 = {2.0};
    std::vector<double> bond_orders = {1.0};
    std::vector<std::vector<double>> hydrogen_data = {atomic_numbers, x_coords, y_coords, z_coords, aid1, aid2, bond_orders};
    Extended_Huckel diatomic_hydrogen(hydrogen_data);
    
    /*
    double hydrogen_enthalpy= diatomic_hydrogen.enthalpy();
    double enthalpy = compound.enthalpy() - hydrogen_enthalpy;
    std::cout << "enthalpy: " << enthalpy << std::endl;
    */

    std::string benzene_filepath = "data/benzene.json";
    std::vector<std::vector<double>>  benzene_data = read_json(benzene_filepath);
    Extended_Huckel benzene(benzene_data);

    double benzene_energy= benzene.total_energy();
    double energy = compound.total_energy() - benzene_energy;
    std::cout << "energy: " << energy << std::endl;
    

    /*
    std::string hydrogen_filepath = "data/ethane.json";
    std::vector<std::vector<double>>  hydrogen_data = read_json(hydrogen_filepath);
    Extended_Huckel diatomic_hydrogen(hydrogen_data);
    double hydrogen_enthalpy = diatomic_hydrogen.enthalpy();
    */

    //std::cout << "H matrix:\n" << compound.H() << std::endl;
    //compound.total_energy();
    //std::cout << "enthalpy: " << compound.enthalpy() << std::endl;
    //std::cout << "Heat summ: " << compound.experimental_heat_summ() << std::endl;
    //std::cout << "E isolation sum: " << compound.E_isol_summ() << std::endl;
    //std::cout << "TOTAL ENERGY:" << std::endl;
    //std::cout << "total energy: " << compound.total_energy() << std::endl;
    
    /*
    // TEST G_AB
    std::cout << compound.G_AB("H", "C", 1.0) << std::endl;
    */

    /*
    // TEST occupation level 
    std::cout << compound.occupation_level(0) << std::endl;
    */

    // TESTS two center off diagonal
    /*
    atomic_orbital u = compound.atomic_orbitals[1];
    atomic_orbital v = compound.atomic_orbitals[4];
    double bohr_radius = 0.52917;
    std::cout<< "u atomic coordinates:\n " << u.coordinates << std::endl;
    std::cout<< "v atomic coordinates:\n " << v.coordinates << std::endl;
    std::cout << "two center off diagonal: " << compound.two_center_off_diagonal(u,v) << std::endl;
    */

    // TESTS for determine_bonda_parameters
    /*
    atomic_orbital u = compound.atomic_orbitals[0];
    atomic_orbital v = compound.atomic_orbitals[9];

    std::cout << "u atom id: " << u.atom_id << std::endl;
    std::cout << "u atom: " << u.atom << std::endl;
    std::cout << "u orbital type: " << u.orbital_type << std::endl;
    std::cout << "u orital axis: " << u.orbital_axis << std::endl;
    std::cout << "u atomic number: " << u.atomic_number << std::endl;
    std::cout << "u coordinates:\n" << u.coordinates << std::endl;

    std::cout << "v atom id: " << v.atom_id << std::endl;
    std::cout << "v atom: " << v.atom << std::endl;
    std::cout << "v orbital type: " << v.orbital_type << std::endl;
    std::cout << "v orital axis: " << v.orbital_axis << std::endl;
    std::cout << "v atomic number: " << v.atomic_number << std::endl;
    std::cout << "v coordinates:\n" << v.coordinates << std::endl;

    std::pair<double, double> bond_paramaters = compound.determine_bond_parameters(u,v,2.0);
    std::cout << "beta: " << bond_paramaters.first << std::endl;
    std::cout << "lambda: " << bond_paramaters.second << std::endl;
    */
    

    
    // TESTS for pp_bond_type
    /*
    atomic_orbital u = compound.atomic_orbitals[1];
    atomic_orbital v = compound.atomic_orbitals[9];
    std::cout << compound.pp_bond_type(u,v) << std::endl;
    */
    

    // TESTS for are_atoms_bonded 
    /*
    atomic_orbital u = compound.atomic_orbitals[5];
    atomic_orbital v = compound.atomic_orbitals[10];
    std::tuple<bool, double> bonded = compound.are_atoms_bonded(u,v);
    std::cout << "bonded: " << std::get<0>(bonded) << std::endl;
    std::cout << "bond order: " << std::get<1>(bonded) << std::endl;
    */

    // TESTS FOR CONTRUCTOR-- atomic_orbitals attribute
    
    /*
    std::vector<atomic_orbital> atomic_orbitals = compound.atomic_orbitals;
    for(const auto& atomic_orbital : atomic_orbitals)
    {
        std::cout << "NEW ORBITAL: " << std::endl;
        std::cout << "atom id: " << atomic_orbital.atom_id << std::endl;
        std::cout << "atom: " << atomic_orbital.atom << std::endl;
        std::cout << "orbital type: " << atomic_orbital.orbital_type << std::endl;
        std::cout << "orital axis: " << atomic_orbital.orbital_axis << std::endl;
        std::cout << "atomic number: " << atomic_orbital.atomic_number << std::endl;
        std::cout << "coordinates:\n" << atomic_orbital.coordinates << std::endl;
    }
    */
    

    // TESTS FOR create_atomic_orbitals
    /*
    arma::vec coordinates = {0.0, 0.0, 0.0};
    std::vector<atomic_orbital> my_orbitals = create_atomic_orbitals(1, "H", 0.0, coordinates, true);
    for (const auto& atomic_orbital : my_orbitals) {
        std::cout << "NEW ORBITAL: " << std::endl;
        std::cout << "atom id: " << atomic_orbital.atom_id << std::endl;
        std::cout << "atom: " << atomic_orbital.atom << std::endl;
        std::cout << "orbital type: " << atomic_orbital.orbital_type << std::endl;
        std::cout << "orital axis: " << atomic_orbital.orbital_axis << std::endl;
        std::cout << "atomic number: " << atomic_orbital.atomic_number << std::endl;
        std::cout << "coordinates:\n" << atomic_orbital.coordinates << std::endl;

    }
    */



   
    return 0;
}




