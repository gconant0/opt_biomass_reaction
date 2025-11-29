//Gavin Conant  12/14/09
//Class for linear programming routinues for flux-balance analysis
#include <iostream>
#include <math.h>
#include <iomanip>
#include <glpk.h> 
#include <gen_dna_funcs.h>
#include "stoich_mat.h"




#ifndef ___LIN_PROGRAM_H___
#define ___LIN_PROGRAM_H___


class Linear_Programmer
{
public:
	Linear_Programmer() {};
	Linear_Programmer(Reaction_matrix *new_mat, int opt_num);
	void set_bounds();
	void do_optimization();
	void do_optimization(bool save_fluxes);
	double get_opt_val();
	void print_soln();
	void save_soln();
	void bound_with_current_opt();
	bool run_success ()	{return(success);};
	Reaction_matrix * get_matrix() {return(curr_matrix);};
	int get_opt_num()  {return(opt_equation_num);};
	~Linear_Programmer() {glp_delete_prob(lp); };
	
protected:
	int non_zero_coeff, *row_indices, *col_indices, opt_equation_num;
	double *coefficients;
	bool success;
	Reaction_matrix *curr_matrix;
	glp_prob *lp; 
	
	
};



#endif













