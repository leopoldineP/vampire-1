#ifndef CREATE_H_
#define CREATE_H_
///
/// @file
/// @brief Contains the cs namespace header. 
///
/// @details This is the detailed description of the funtion of this file
///
/// @section notes Implementation Notes
/// This is a list of other notes, not related to functionality but rather to implementation. 
/// Also include references, formulae and other notes here.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
#include <vector>
#include <cmath>

/// @namespace
/// @brief Contains all functions and data associated with system creation in zspin.
/// 
/// @internal
///=====================================================================================
///

namespace cs{
	class catom_t {
		public:
			
			// Coordinates
			double x;
			double y;
			double z;
		
			// Flags
			bool include;

			// Integers
			int material;
			int uc_category;
			int lh_category;
			int grain;
			int supercell;
			int mpi_type;
			int mpi_cpuid;
			int mpi_atom_number;
	  int mpi_old_atom_number;
			catom_t():
				x(0.0),
				y(0.0),
				z(0.0),
				include(false),
				material(0),
				uc_category(0),
				lh_category(0),
				grain(0),
				supercell(0),
				mpi_type(0),
				mpi_cpuid(0),
				mpi_atom_number(0),
				mpi_old_atom_number(0)
			{};
};
/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int create();

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int create_crystal_structure(std::vector<cs::catom_t> &);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int create_system_type(std::vector<cs::catom_t> &);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int create_neighbourlist(std::vector<cs::catom_t> &, std::vector<std::vector <int> > &);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int set_atom_vars(std::vector<cs::catom_t> &, std::vector<std::vector <int> > &);

int voronoi_film(std::vector<cs::catom_t> &);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int bulk(std::vector<cs::catom_t> &,const int);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int cube(double[], std::vector<cs::catom_t> &,const int);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int sphere(double[], std::vector<cs::catom_t> &,const int);
	

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int cylinder(double[], std::vector<cs::catom_t> &,const int);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int truncated_octahedron(double[], std::vector<cs::catom_t> &,const int);

int sort_atoms_by_grain(std::vector<cs::catom_t> &);
int clear_atoms(std::vector<cs::catom_t> &);
}

#endif /*CREATE_H_*/
