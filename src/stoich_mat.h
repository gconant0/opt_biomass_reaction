#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include "gen_dna_funcs.h"
#include <list>
#include "liblapack.h"


#ifndef ___STOICH_MAT_H___
#define ___STOICH_MAT_H___



enum BOUND_TYPE {BD_FREE, BD_UPPER, BD_LOWER, BD_BOTH, BD_FIXED};

enum RXN_RELATE {PLANE, LINE, EMPTY};
class Metabolite 
{
public:
	int met_num;
	std::string met_name;
	bool visit;
	int operator==(Metabolite test);
	Metabolite& operator= (Metabolite &assign_from);
};

class Reaction
{
public:
	int react_num, num_react, num_prod;
	double *react_stoich, *prod_stoich, upper_bound, lower_bound, old_lower, old_upper, flux, min, max, base_flux;
    std::string react_name;
	Metabolite **reactants, **products;
	BOUND_TYPE type, old_type;
	bool reversible, exchange, lethal;
	
	Reaction();
	Reaction& operator= (Reaction& assign_from);
	void set_products_reactants(std::list<Metabolite*> &reactant_list, std::list<Metabolite*> &product_list, std::list<double> &r_stoich, std::list<double> &p_stoich);
    void knockout_rxn();
    void restore_rxn();
    ~Reaction();
};



class Reaction_matrix 
{
public:
	int dim_null_space;
	int *nrhs, *fortran_cols, *fortran_rows, *info, *lwork, *pivots, *iwork;
	char *trans;
	double *fortran_matrix, *sing_vals, *u_mat, *v_mat,  *work,
	**null_space_span;
   
	
	Reaction_matrix();
	Reaction_matrix(std::string infile);
    Reaction_matrix(std::istream& datass);
	Reaction_matrix& operator= (Reaction_matrix& assign_from);
	int get_num_reactions()         {return(num_reactions);};
	int get_num_metabolites()		{return(num_metabolites);};
	void create_fortran_stoic_matrix();
	void compute_ortho_basis();
	Metabolite* get_metabolite(int i)  {return(metabolites[i]);};
    Metabolite* get_metabolite_by_name(std::string name);
    Reaction* get_reaction_by_name(std::string name);
	Reaction* get_reaction(int i)   {return(reactions[i]);};
	int get_omit_reaction()			{return(omit_reaction);};
	void set_omit_reaction(int s);
	~Reaction_matrix();
	void real_omit_reaction(int s);
	
protected:
	int num_reactions, num_metabolites, omit_reaction, omit_size, knockout_rxn;
	bool internal_data;
	Metabolite **metabolites;
	Reaction **reactions;
    
    std::map<std::string, Metabolite*> Metabolite_hash;
    std::map<std::string, Reaction*> Reaction_hash;
    void read_reaction_matrix(std::istream& fin);
};


class Space_compare {
public:
	int cnt_preclude;
	int *fortran_cols, *fortran_rows, *nrhs, *info, *lwork;
	char *trans;
	double *fortran_matrix, *fortran_vecs, *work;
	
	Space_compare();
	double do_comparison (Reaction_matrix *full_matrix, Reaction_matrix *reduced_matrix);
	~Space_compare();
};

RXN_RELATE paired_rxn_line_soln(Reaction_matrix * the_matrix, int rxn1, int rnx2);

void get_bounds (Reaction_matrix *curr_matrix, std::string bound_file);
void get_bounds (Reaction_matrix *curr_matrix, std::stringstream& bound_ss);
void get_bounds_read (Reaction_matrix *curr_matrix, std::istream& fin);
void get_bounds_var_search (Reaction_matrix *curr_matrix, std::string bound_file, bool *&b_variable);

#endif
