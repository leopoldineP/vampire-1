//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package ... ?
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <fstream>
#include <algorithm>

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
    bool y = false;
    bool z = false;
};

// Function forward declaration
std::vector < seed_point_defects > generate_random_defect_seed_points(int defect_amount);
void voronoi_defects(std::vector<cs::catom_t> & catom_array);


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

    // Print informative message to logfile
    zlog << zTs() << "Calculating defect properties of system" << std::endl;


    //decide which shape 
    bool sphere = false; //hardcoded as true for now, change to a different version in interface 
    bool triangle = false;
    if (sphere==true){ 
        // Print informative message to screen
        zlog << zTs() << "Creating spherical defects" << std::endl;

        //call function to calculate positions of defects 
        std::vector < seed_point_defects > defect_pos = generate_random_defect_seed_points(defect_amount);

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
                    zlog << zTs() << "delete" << atom << std::endl;
                }
            }
        }
    }//end of sphere

    else if (triangle==true){ //needs new algorithm
        for (int def=0;def<create::internal::defect_amount;def++){

            // Print informative message to screen
            zlog << zTs() << "Creating triangular prism defects" << std::endl;

            //call function to calculate positions of defects 
            std::vector < seed_point_defects > defect_pos = generate_random_defect_seed_points(defect_amount);

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
         voronoi_defects(catom_array);
    
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

    //open file for distribution testing
    std::ofstream distributionfile;
    distributionfile.open ("distribution.txt");

    //set space for random number generator
    const double defectspace_x = abs(defectspace_max_x)-abs(defectspace_min_x);
	const double defectspace_y = abs(defectspace_max_y)-abs(defectspace_min_y);
    const double defectspace_z = abs(defectspace_max_z)-abs(defectspace_min_z);

    //set distribution of random points
    defect_distribution reversed;
    defect_distribution gaussian;
    //defect_distribution uniform; - default setting
    
    //parameters - move to user interface
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
            position.x = (mtrandom::gaussian())*defectspace_x+defectspace_min_x;
        }
	    else{ //uniform
            position.x = (create::internal::grnd())*defectspace_x+defectspace_min_x;
        }
        //y-coordinate
        if (gaussian.y==true){
            position.y = (mtrandom::gaussian())*defectspace_y+defectspace_min_y;
        }
	    else{ //uniform
            position.y = (create::internal::grnd())*defectspace_y+defectspace_min_y;
        }
        //z-coordinate
        if (gaussian.z==true){
            position.z = (mtrandom::gaussian())*defectspace_z+defectspace_min_z;
        }
	    else{ //uniform
            position.z = (create::internal::grnd())*defectspace_z+defectspace_min_z;
        }

        //reversed gaussian distribution
        if (reversed.x==true){
            if (position.x-defectspace_min_x<defectspace_x/2){
                //if in first half, shift to second half, this splits the gaussian curve into two and shifts its center away from the sample's center to the edges 
                position.x=position.x+defectspace_x/2;
            }
            if (position.x-defectspace_min_x>defectspace_x/2){
                //if in second half, shift to first half 
                position.x=position.x-defectspace_x/2;
            }
        }
        if (reversed.y==true){
            if (position.y-defectspace_min_y<defectspace_y/2){
                //if in first half, shift to second half 
                position.y=position.y+defectspace_y/2;
            }
            if (position.y-defectspace_min_y>defectspace_y/2){
                //if in second half, shift to first half 
                position.y=position.y-defectspace_y/2;
            }
        }
        if (reversed.z==true){
            if (position.z-defectspace_min_z<defectspace_z/2){
                //if in first half, shift to second half 
                position.z=position.z+defectspace_z/2;
            }
            if (position.z-defectspace_min_z>defectspace_z/2){
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

        //distribution check: write defects to file 
        distributionfile << position.x << " " << position.y << " " << position.z << ".\n";
        
    }
    distributionfile << "space " << defectspace_x << " " << defectspace_y << " " << defectspace_z << ".\n";
    
    

    distributionfile.close();

    return defects_positions;

} //end of defects position function


void voronoi_defects(std::vector<cs::catom_t> & catom_array){
    //--------------------------------------------------------------------------------------------
    // part 1: create voronoi structure in chosen area
    //--------------------------------------------------------------------------------------------

    //make a copy of catom_array with only relevant atoms (identical to catom_array if defectspace=systemsize so for now not included) 
    
    uint64_t num_local_atoms = catom_array.size();
    uint64_t num_total_atoms = 0; // number of atoms across all processors

      #ifdef MPICF
         // calculate number of atoms to be output on all processors
         MPI_Allreduce(&num_local_atoms, &num_total_atoms, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
      #else
         num_total_atoms = num_local_atoms;
      #endif
    int nc = floor(num_total_atoms/(defect_vacancies+0.5*defect_vacancies)); //number of voronoi cells to be generated in space (currently have on average 50% more atoms in cell than vacancies needed)

    std::vector< seed_point_defects > voronoiseed_positions(0);//vector to store all positions of voronoi seeds

    //loop through cell amount to create random position for each (uniformly distributed)
    for (int i=0; i<nc; i++){

        // generate random x,y,z trial point 
		seed_point_defects voronoiseed; 
        
        voronoiseed.x = (create::internal::grnd()*cs::system_dimensions[0]); //x-coordinate 
        voronoiseed.y = (create::internal::grnd()*cs::system_dimensions[1]); //y-coordinate
        voronoiseed.z = (create::internal::grnd()*cs::system_dimensions[2]); //z-coordinate

        //test positions for being good choices?

        //store positions
        voronoiseed_positions.push_back(voronoiseed);
    }

    std::vector< std::vector <int> > voronoigrid(nc-1); //vector storing atoms that correspond to a cell
    std::vector< std::vector <double> > sortingdistances(nc-1); //temporary vector storing distances of the atom to their seed point of the cell for sorting
    //loop through atoms inside defectspace to assign to a seed point they're closest to
    for (int atom=0;atom<catom_array.size();atom++){

        //find correct voronoi cell
        std::vector< double > seeddistances(0); 
        for (int i=0;i<nc;i++){
            //calculate distance of atom from the seed position
            double distance_from_seed_sq= (catom_array[atom].x-voronoiseed_positions[i].x)*(catom_array[atom].x-voronoiseed_positions[i].x) + 
                                          (catom_array[atom].y-voronoiseed_positions[i].y)*(catom_array[atom].y-voronoiseed_positions[i].y) + 
                                          (catom_array[atom].z-voronoiseed_positions[i].z)*(catom_array[atom].z-voronoiseed_positions[i].z);
            //save all distances temporarily
            seeddistances.push_back(distance_from_seed_sq);
        }
        //find closest seed point
        int closestseedindex = std::min_element(seeddistances.begin(),seeddistances.end()) - seeddistances.begin();
        
        //add to that cell
        voronoigrid[closestseedindex].push_back(atom);
        sortingdistances[closestseedindex].push_back(seeddistances[closestseedindex]);
    }

    //sort voronoigrid vector (bubble sort)
    for (int cells=0;cells<nc;cells++){
        for (int j=0;j<sortingdistances[cells].size();j++){
            for (int i=j+1;i<sortingdistances[cells].size();i++){
                if(sortingdistances[cells][j] < sortingdistances[cells][i]){
                    //sort atoms accordingly
                    int temp = voronoigrid[cells][i];
                    voronoigrid[cells][i] = voronoigrid[cells][j];
                    voronoigrid[cells][j] = temp;
                    //also update distances
                    int temp_d = sortingdistances[cells][i];
                    sortingdistances[cells][i] = sortingdistances[cells][j];
                    sortingdistances[cells][j] = temp_d;
                }
            }
        }
    }
    //end of distances vector

    //check grid - comment out or delete
    if (vmpi::my_rank==0){
        std::ofstream voronoigridfile;
    voronoigridfile.open ("voronoigrid.txt");
    //head of file: cell, amount of atoms in cell
    for (int h=0;h<nc;h++){
        voronoigridfile << h << " " << voronoigrid[h].size() << ".\n";
    }
    //rest of file filled with voronoicell, atomindex, atom coordinates xyz
    for (int h=0;h<nc;h++){
        for (int k=0;k<voronoigrid[h].size();k++){
            voronoigridfile << h << " " << voronoigrid[h][k] << " " << catom_array[voronoigrid[h][k]].x << " " << catom_array[voronoigrid[h][k]].y << " " << catom_array[voronoigrid[h][k]].z << ".\n";
        }
    }
    voronoigridfile.close();
    //end of checks
    }


    //--------------------------------------------------------------------------------------------
    // part 2: create defects from structure
    //--------------------------------------------------------------------------------------------
    int loopbreak=0; //variable to break the loop from running if a certain number of failed trial positions has been reached
    std::vector< int > chosencells(0);//temporary vector for checking

    //open file for distribution testing
    std::ofstream voronoidistributionfile;
    voronoidistributionfile.open ("voronoidistribution.txt");

    for (int def=0;def<defect_amount;def++){
        int defectseed = floor(create::internal::grnd()*nc);
        //mpi_broadcast defectseed;
        zlog << zTs() << "defectseed " << defectseed << std::endl;
        
        //check that this voronoi cell has not been chosen as a defect before
        bool newcell=true;
        for (int i=0;i<chosencells.size();i++){
            if (chosencells[i]==defectseed || voronoigrid[defectseed].size()<defect_vacancies){
                newcell=false;
                break;
            }
        } //check complete

        if (newcell){ //delete a portion of/all atoms in this voronoi cell
            chosencells.push_back(defectseed);
            for(int vac=0;vac<defect_vacancies;vac++){ 
                if (vac<voronoigrid[defectseed].size()){
                    catom_array[voronoigrid[defectseed][vac]].include=false; //delete atom from catom array
                    //work out how many atoms are being deleted on each processor 
                    //you can bring all atom positions and radii from center, source processor number + radius
                    // mpi gather 
                    //std::vector <int> local_proc(sortingdistances.size() of current defect,vmpi::my_rank);
                    //1. reduction: how many atoms in this grain =num_atoms_grain
                    //std::vector <int> processor_of_atom(num_atoms_grain);
                    //vmpi::collate(local_proc, processor_of_atom)
                    //print: should be list of total atoms in defect and where each one is
                    //collate on radii(local_radii, global_radii)
                    //sort by radius, iterate until vacancies are reached 

                    //distribution check: write defects to file 
                    voronoidistributionfile << def << " " << vac << " " << catom_array[voronoigrid[defectseed][vac]].x << " " << catom_array[voronoigrid[defectseed][vac]].y << " " << catom_array[voronoigrid[defectseed][vac]].z << ".\n";
                } 
                else { 
                    zlog << zTs() << "Maximum amount of atoms in voronoi cell was reached at " << vac << " for defect " << def << std::endl;
                    break; //continue with next defect -alternative to be coded: take defects from neighbouring cells
                } 
            }
            loopbreak=0; //reset loopbreak
        }
        else { //cell chosen was not valid
            def=def-1;
            loopbreak=loopbreak+1;
            if (loopbreak>1000000){
                // Print informative message to log book
                zlog << zTs() << "Maximum amount of trial positions was reached for defect " << def+1 << " and simulation was ended" << std::endl;
                err::vexit(); //end simulation
            }
        }
    }

    voronoidistributionfile.close();


} //end of voronoi function

} // end of internal namespace
} // end of create namespace