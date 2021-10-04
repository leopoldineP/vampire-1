//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package ... ?
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <fstream>

// Vampire headers
#include "create.hpp"

// Internal create header
#include "internal.hpp"
#include "vio.hpp"

namespace create{

namespace internal{


//define structure to store positions of defects
struct seed_point_defects{
    double x; //x-position
    double y; //y-position
    double z; //z-position
};

//define structure to store area parameters for enclosed area in which defects can be found
struct area{
    double min_x; //start of defined area in x
    double max_x; //end of defined area in x
    double min_y; //start of defined area in y
    double max_y; //end of defined area in y
};

// Function forward declaration
std::vector < seed_point_defects > generate_random_defect_seed_points(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int defect_amount);


//-----------------------------------------------------------------------------
//
// Function to create defects in specified shape
// around previously defined sites of defects 
//
// (c) 
//
//-----------------------------------------------------------------------------
void defects (std::vector<cs::catom_t> & catom_array){

    // Defect properties not guaranteed to exist
	// return here if unused to avoid segmentation fault
	//if(create::internal::mp.size() == 0) return; -copied from alloy, need to change?

    // Print informative message to screen
    zlog << zTs() << "Calculating defect properties of system" << std::endl;

    //parameteres (- for user interface compatability mp::num_materials;)
    const int defect_amount = 5; //local constant for number of defects 
    const int defect_vacancies = 10; //how many sites are missing per defect


    //decide which shape 

    //if (sphere=true){ 
      //  for (){

       // Print informative message to screen

            //delete points inside shape (extension: only delete a certain percentage of points at the edge)
        //    catom_array[atom].include=false;

         
        //}
       
    //}

     //else { //irregular shape (default version)

      // Print informative message to screen
    

     //}

     return;

} //end of defect shape function


//-----------------------------------------------------------------------------
//
// Function to generate random positions for defects 
// or read in manuallly set positions from .dcf file
//
// (c)
//
//-----------------------------------------------------------------------------
std::vector < seed_point_defects > generate_random_defect_seed_points(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int defect_amount){

    //parameters
    double min_distance = 0; //minimum distance required between defects - none by default / defects should not be touching ?
    double max_distance = sizex+sizey+sizez; //maximum distance required between defects - none by default
    
    //set space if defects are in a defined area - by default equal to system dimensions
    const double sizex = cs::system_dimensions[0];
	const double sizey = cs::system_dimensions[1];
    const double sizez = cs::system_dimensions[2];
    //read in sizes from input file for defined space (area struct)

    //vector for defect positions
    std::vector< seed_point_defects > defects(0);

    // re-seed random number generator on each CPU 
	create::internal::grnd.seed(vmpi::parallel_rng_seed(create::internal::defect_seed));

    //loop through defects to create random position for each
    for (int i=0; i<defect_amount; i++){

        // generate random x,y,z trial point
		seed_point_defects position;
	    position.x = (create::internal::grnd()*1.4-0.2)*size_x;
	    position.y = (create::internal::grnd()*1.4-0.2)*size_y;
	    position.z = (create::internal::grnd()*1.4-0.2)*size_z;

		// flag to see if positions are too close to each other (and fullfill specified min/max restrictions)
		bool defect_distance=false;

		// loop over all previous positions and check if position is valid within restrictions placed on the position relative to other positions
		for(unsigned int g=0; g<defects.size(); g++){
			double distance_x = position.x-defects[g].x;
			double distance_y = position.y-defects[g].y;
            double distance_z = position.z-defects[g].z;
			double distance_ij = sqrt(distance_x*distance_x + distance_y*distance_y + distance_z*distance_z);
			if(distance_ij<min_distance || distance_ij>max_distance){ //replace with && if minimum AND maximum are chosen, rn either or
				defect_distance = true;
				break;
			}
		}

		// save valid positions
		if(defect_distance == false) defects.push_back(position);
    }

    return defects;

} //end of defects position function

} // end of internal namespace
} // end of create namespace