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


namespace create{

namespace internal{


//define structure to store positions of defects
struct seed_point_defects{
    double x; //x-position
    double y; //y-position
    double z; //z-position
}

//define structure to store area parameters for enclosed area in which defects can be found
struct area{
    double min_x; //start of defined area in x
    double max_x; //end of defined area in x
    double min_y; //start of defined area in y
    double max_y; //end of defined area in y
}

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


// //end of defect shape function


//-----------------------------------------------------------------------------
//
// Function to generate random positions for defects 
// or read in manuallly set positions from .dcf file
//
// (c)
//
//-----------------------------------------------------------------------------
std::vector < seed_point_defects > generate_random_defect_seed_points(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int defect_amount){
    
    //set space if defects are in a defined area

    //vector for defect positions
    std::vector< seed_point_defects > defects(0);

    
    return defects;

} //end of defects position function

} // end of internal namespace
} // end of create namespace