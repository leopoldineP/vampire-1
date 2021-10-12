//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package ... ?
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <fstream>

// Vampire headers
#include "create.hpp"
#include "errors.hpp"

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
std::vector < seed_point_defects > generate_random_defect_seed_points(const int defect_amount);


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
    const int defect_amount = 2; //local constant for number of defects 

    //call function to calculate positions of defects
    std::vector < seed_point_defects > defect_pos = generate_random_defect_seed_points(defect_amount);


    //decide which shape 

    //if (sphere=true){ 
    for (int def=0;def<defect_amount;def++){

        // Print informative message to screen
        zlog << zTs() << "Creating spherical defects" << std::endl;

        //set defects radius
        double defect_radius_squared = 2.0;

        //loop over all atoms to see what atoms are within sphere - is there a quicker way to do this? (neighbourlist of position?)
        int number_of_atoms = catom_array.size();
        for (int atom=0;atom<number_of_atoms;atom++){

            //calculate distance of atom from the defect position
            double distance_from_defect_sq= (catom_array[atom].x-defect_pos[def].x)*(catom_array[atom].x-defect_pos[def].x) + 
                                            (catom_array[atom].y-defect_pos[def].y)*(catom_array[atom].y-defect_pos[def].y) + 
                                            (catom_array[atom].z-defect_pos[def].z)*(catom_array[atom].z-defect_pos[def].z);

            if (distance_from_defect_sq<=defect_radius_squared){
                //delete points inside shape (extension: only delete a certain percentage of points at the edge)
                catom_array[atom].include=false;
            }
        }
    }
       
    //}

     //else { //irregular shape (default version)

      // Print informative message to screen
    

     //}

     //save defect attributes to file ?

     

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
std::vector < seed_point_defects > generate_random_defect_seed_points(const int defect_amount){

    //set space if defects are in a defined area - by default equal to system dimensions
    const double defectspace_x = cs::system_dimensions[0];
	const double defectspace_y = cs::system_dimensions[1];
    const double defectspace_z = cs::system_dimensions[2];
    //read in sizes from input file for defined space (area struct)

    //parameters - move to user interface
    double min_distance = 0; //minimum distance required between defects - none by default / defects should not be touching ?
    double max_distance = defectspace_x + defectspace_y + defectspace_z; //maximum distance required between defects - none by default
    int max_trial_positions = 1000;
    
    //vector for defect positions
    std::vector< seed_point_defects > defects_positions(0);

    // re-seed random number generator on each CPU 
	create::internal::grnd.seed(vmpi::parallel_rng_seed(create::internal::defect_seed));

    //loop through defects to create random position for each
    for (int i=0; i<defect_amount; i++){

        // generate random x,y,z trial point
		seed_point_defects position;
	    position.x = (create::internal::grnd()*1.4-0.2)*defectspace_x;
	    position.y = (create::internal::grnd()*1.4-0.2)*defectspace_y;
	    position.z = (create::internal::grnd()*1.4-0.2)*defectspace_z;

		// flag to see if positions are too close to each other (and fullfill specified min/max restrictions)
		bool defect_distance=true;

		// loop over all previous positions and check if position is valid within restrictions placed on the position relative to other positions
        int check_loop=0; //counter to stop program from being stuck in this loop
		for (unsigned int g=0; g<defects_positions.size(); g++){
			double distance_x = position.x-defects_positions[g].x;
			double distance_y = position.y-defects_positions[g].y;
            double distance_z = position.z-defects_positions[g].z;
			double distance_ij = sqrt(distance_x*distance_x + distance_y*distance_y + distance_z*distance_z);
			if(distance_ij<min_distance || distance_ij>max_distance){ //replace with && if minimum AND maximum are chosen, rn either or
				defect_distance = false;
                i=i-1; //run the random seed generator again for this defect
                check_loop=check_loop+1; //counting unsuccessful trial positions
                break;
			}
		}

        //stop  if maximum of unsuccessful trial positions has been reached 
        if (check_loop>max_trial_positions){
            // Print informative message to screen
            zlog << zTs() << "Maximum amount of trial positions was reached for defect " << i+1 << "and simulation was ended" << std::endl;
            err::vexit(); //end simulation
        }

		// save valid positions
		if(defect_distance == true) defects_positions.push_back(position);
       
    }

    return defects_positions;

} //end of defects position function

} // end of internal namespace
} // end of create namespace