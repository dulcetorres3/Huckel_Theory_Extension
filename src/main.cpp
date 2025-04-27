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
    int atom_id;
    std::string atom;
    std::string orbital_type;
    double atomic_number;
    arma::vec coordinates;

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

class Extended_Huckel
{
    private:
    int N;
    //std::vector<atomic_orbital> atomic_orbitals;
    
    public:
    std::vector<atomic_orbital> atomic_orbitals;
    Extended_Huckel(std::vector<std::vector<double>> compound_data)
    {
        std::vector<double> atomic_numbers = compound_data[0];
        std::vector<double> x_coords = compound_data[1];
        std::vector<double> y_coords = compound_data[2];
        std::vector<double> z_coords = compound_data[3];

        int total_atoms = atomic_numbers.size();
        for(int i=0; i < total_atoms; i++)
        {
            double atomic_number = atomic_numbers[i];
            std::string atom = element_map[atomic_number];
            int atom_id = i;
            arma::vec coordinates = {x_coords[i], y_coords[i], z_coords[i]};

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

    /*
    double two_center_off_diagonal(atomic_orbital u, atomic_orbital v)
    {

    }
    arma::mat H_matrix()
    {
        arma::mat H(N,N);

        for(int u=0; u < N; u++)
        {
            atomic_orbital u_orbital = atomic_orbitals[u];
            int u_atom_id = u_orbital.atom_id;
            for(int v=0; v < N; v++)
            {
                double uv;
                atomic_orbital v_orbital = atomic_orbitals[v];
                int v_atom_id = v_orbital.atom_id;

                // diagonal terms
                if(u==v)
                {

                }
                else
                {
                    // one-center off diagonal terms
                    if(u_atom_id==v_atom_id)
                    {
                        uv=0.0;
                    }
                    else // two-center off diagonal terms
                    {

                    }
                }

            }
        }


    }
*/

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
    
    // Extract elements and coordinates
    std::vector<double> elements;
    std::vector<double> x_coords, y_coords, z_coords;
    
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

    std::vector<std::vector<double>>  compound_data = {elements, x_coords, y_coords, z_coords};
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




