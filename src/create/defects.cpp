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
#include "random.hpp"

namespace create{

namespace internal{


//define structure to store positions of defects
struct seed_point_defects{
    double x; //x-position
    double y; //y-position
    double z; //z-position
};

//define structure for distribution checks
struct defect_distribution{
    bool x = false;
    bool y = true;
    bool z = true;
};

//define structure to store area parameters for enclosed area in which defects can be found
struct area{
    double min_x; //start of defined area in x
    double max_x; //end of defined area in x
    double min_y; //start of defined area in y
    double max_y; //end of defined area in y
    double min_z; //start of defined area in z
    double max_z; //end of defined area in z
};

// Function forward declaration
std::vector < seed_point_defects > generate_random_defect_seed_points(int defect_amount);
void voronoi_defects(std::vector<cs::catom_t> & catom_array, int defect_amount);


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

    //call function to calculate positions of defects -has to be inside function as not valid for voronoi
    std::vector < seed_point_defects > defect_pos = generate_random_defect_seed_points(defect_amount);


    //decide which shape 

    bool sphere = true; //hardcoded as true for now, change to a different version in interface 
    bool triangle = false;

    if (sphere==true){ 
        // Print informative message to screen
        zlog << zTs() << "Creating spherical defects" << std::endl;

        for (int def=0;def<create::internal::defect_amount;def++){ //loop through defects

            //set defects radius
            double defect_radius_squared = 40.0;
            
            //add defects plus position to log file 
            zlog << zTs() << "d " << def << " " << defect_pos[def].x << " " << defect_pos[def].y << " " << defect_pos[def].z << std::endl;


            //loop over all atoms to see what atoms are within sphere - is there a quicker way to do this? (neighbourlist of position?)
            int number_of_atoms = catom_array.size();
            for (int atom=0;atom<number_of_atoms;atom++){

                //calculate distance of atom from the defect position
                double distance_from_defect_sq= (catom_array[atom].x-defect_pos[def].x)*(catom_array[atom].x-defect_pos[def].x) + 
                                                (catom_array[atom].y-defect_pos[def].y)*(catom_array[atom].y-defect_pos[def].y) + 
                                                (catom_array[atom].z-defect_pos[def].z)*(catom_array[atom].z-defect_pos[def].z);

                if (distance_from_defect_sq<=defect_radius_squared){
                    //make atoms inside shape nonmagnetic (extension: only delete a certain percentage of atoms in the centre)
                    catom_array[atom].include=false; //delete
                    //mp::material[catom_array[atom].material].non_magnetic = 1; //make non-magnetic
                    zlog << zTs() << "d" << atom << std::endl;
                }
            }
        }
    }//end of sphere

    else if (triangle==true){ //needs new algorithm
        for (int def=0;def<create::internal::defect_amount;def++){

            // Print informative message to screen
            zlog << zTs() << "Creating triangular prism defects" << std::endl;

            //set defects size
            double deftriangle_height = 2.0; //hard coded but take from interface in the future
            double deftriangle_base = 2.0; //change
            double deftriangle_area =(deftriangle_height*deftriangle_base)/2.0; 
            bool insidetriangle = false; //use barycentric coordinates instead
       
            //loop over all atoms to see what atoms are within triangular prism - is there a quicker way to do this? (neighbourlist of position?)
            int number_of_atoms = catom_array.size();
            for (int atom=0;atom<number_of_atoms;atom++){

                //calculate corners of triangle
                

                if (insidetriangle == true){
                    //make atoms inside shape nonmagnetic (extension: only delete a certain percentage of atoms in the centre)
                    // delete: catom_array[atom].include=false;
                    mp::material[catom_array[atom].material].non_magnetic = 1;
                }
            }
        }
    }//end of triangle

    else { //irregular shape:voronoi (default version)

         // Print informative message to screen
         zlog << zTs() << "Creating voronoi defects" << std::endl;

         //call voronoi function
         voronoi_defects(catom_array,defect_amount);
    

     } //end of irregular voronoi defects

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
std::vector < seed_point_defects > generate_random_defect_seed_points(int defect_amount){

    //read in sizes from input file for defined space - by default equal to system dimensions
    area defined_space;
    defined_space.min_x = 0;
    defined_space.max_x = cs::system_dimensions[0];
    defined_space.min_y = 0;
    defined_space.max_y = cs::system_dimensions[1];
    defined_space.min_z = 0;
    defined_space.max_z = cs::system_dimensions[2];

    //set space for random number generator
    const double defectspace_x = abs(defined_space.max_x)-abs(defined_space.min_x);
	const double defectspace_y = abs(defined_space.max_y)-abs(defined_space.min_y);
    const double defectspace_z = abs(defined_space.max_z)-abs(defined_space.min_z);

    //set distribution of random points
    defect_distribution reversed;
    defect_distribution gaussian;
    //defect_distribution uniform; - default setting
    
    //parameters - move to user interface
    //minimum distance required between defects - none by default / defects should not be touching ?
    double max_distance = defectspace_x + defectspace_y + defectspace_z; //maximum distance required between defects - none by default
    int max_trial_positions = 1000;
    
    //vector for defect positions
    std::vector< seed_point_defects > defects_positions(0);

    // re-seed random number generator on each CPU 
	create::internal::grnd.seed(vmpi::parallel_rng_seed(create::internal::defect_seed));

    //loop through defects to create random position for each
    for (int i=0; i<create::internal::defect_amount; i++){

        // generate random x,y,z trial point 
		seed_point_defects position;
        //x-coordinate
        if (gaussian.x==true){
            position.x = (mtrandom::gaussian())*defectspace_x+defined_space.min_x;
        }
	    else{ //uniform
            position.x = (create::internal::grnd()*1.4-0.2)*defectspace_x+defined_space.min_x;
        }
        //y-coordinate
        if (gaussian.y==true){
            position.y = (mtrandom::gaussian())*defectspace_y+defined_space.min_y;
        }
	    else{ //uniform
            position.y = (create::internal::grnd()*1.4-0.2)*defectspace_y+defined_space.min_y;
        }
        //z-coordinate
        if (gaussian.z==true){
            position.z = (mtrandom::gaussian())*defectspace_z+defined_space.min_z;
        }
	    else{ //uniform
            position.z = (create::internal::grnd()*1.4-0.2)*defectspace_z+defined_space.min_z;
        }

        //reversed gaussian distribution
        if (reversed.x==true){
            if (position.x-defined_space.min_x<defectspace_x/2){
                //if in first half, shift to second half, this splits the gaussian curve into two and shifts its center away from the sample's center to the edges 
                position.x=position.x+defectspace_x/2;
            }
            if (position.x-defined_space.min_x>defectspace_x/2){
                //if in second half, shift to first half 
                position.x=position.x-defectspace_x/2;
            }
        }
        if (reversed.y==true){
            if (position.y-defined_space.min_y<defectspace_y/2){
                //if in first half, shift to second half 
                position.y=position.y+defectspace_y/2;
            }
            if (position.y-defined_space.min_y>defectspace_y/2){
                //if in second half, shift to first half 
                position.y=position.y-defectspace_y/2;
            }
        }
        if (reversed.z==true){
            if (position.z-defined_space.min_z<defectspace_z/2){
                //if in first half, shift to second half 
                position.z=position.z+defectspace_z/2;
            }
            if (position.z-defined_space.min_z>defectspace_z/2){
                //if in second half, shift to first half 
                position.z=position.z-defectspace_z/2;
            }
        }

	
		// flag to see if positions are too close to each other (and fullfill specified min/max restrictions)
		bool defect_distance=true;

		// loop over all previous positions and check if position is valid within restrictions placed on the position relative to other positions
        int check_loop=0; //counter to stop program from being stuck in this loop
		for (unsigned int g=0; g<defects_positions.size(); g++){
			double distance_x = position.x-defects_positions[g].x;
			double distance_y = position.y-defects_positions[g].y;
            double distance_z = position.z-defects_positions[g].z;
			double distance_ij = sqrt(distance_x*distance_x + distance_y*distance_y + distance_z*distance_z);
			if(distance_ij<create::internal::min_defect_distance || distance_ij>max_distance){ //replace with && if minimum AND maximum are chosen, rn either or
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


void voronoi_defects(std::vector<cs::catom_t> & catom_array, int defect_amount){

    

} //end of voronoi function

} // end of internal namespace
} // end of create namespace