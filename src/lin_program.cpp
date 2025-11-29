
#include "lin_program.h"

using namespace::std;


Linear_Programmer::Linear_Programmer(Reaction_matrix *new_mat, int opt_num)
{
	int i, j, cnt;
	char row_name[7], col_name[7], dummy[100];
	
	opt_equation_num=opt_num;
	
	curr_matrix=new_mat;
	
	
	lp = glp_create_prob(); 
	glp_set_prob_name(lp, "FBA"); 
	glp_set_obj_dir(lp, GLP_MAX); 
	
	//Create a row for each metabolite
	glp_add_rows(lp, curr_matrix->get_num_metabolites()); 
	
	//Set the constraints for those rows to fixed 0
	for(i=1; i<=curr_matrix->get_num_metabolites(); i++) {
		strcpy(row_name, "M");
		int_to_string(dummy, 99, i);
		strcat(row_name, dummy);
       // cout<<"Row "<<i<<": "<<row_name<<endl;
		glp_set_row_name(lp, i, row_name); 
		glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
	}
	
	
	//Create a column for each reaction
	glp_add_cols(lp, curr_matrix->get_num_reactions()); 
	
	for(i=1; i<=curr_matrix->get_num_reactions(); i++) {
		strcpy(col_name, "R");
		int_to_string(dummy, 99, i);
		strcat(col_name, dummy);
		
        //cout<<"Column "<<i<<": "<<col_name<<endl;
        
		glp_set_col_name(lp, i, col_name);	
		
		if (i == (opt_equation_num+1))
			glp_set_obj_coef(lp, i, 1.0);
		else 
			glp_set_obj_coef(lp, i, 0.0);
	}
	set_bounds();
		
	non_zero_coeff=0;
	
    cout<<"NM: "<<curr_matrix->get_num_metabolites()<<" NR: "<<curr_matrix->get_num_reactions()<<endl;
	
	for(j=0; j<curr_matrix->get_num_reactions(); j++) {
		for(i=0; i<curr_matrix->get_reaction(j)->num_react; i++) 
			non_zero_coeff++;
		for(i=0; i<curr_matrix->get_reaction(j)->num_prod; i++) 
			non_zero_coeff++;
	}
			
	
	row_indices=new int [non_zero_coeff+1];
	col_indices=new int [non_zero_coeff+1];
	coefficients=new double [non_zero_coeff+1];
	cnt=0;
	
	for(j=0; j<curr_matrix->get_num_reactions(); j++) {
		for(i=0; i<curr_matrix->get_reaction(j)->num_react; i++) {
			cnt++;
			row_indices[cnt]=curr_matrix->get_reaction(j)->reactants[i]->met_num+1;
			col_indices[cnt]=j+1;
			coefficients[cnt]=(double)(curr_matrix->get_reaction(j)->react_stoich[i]);
		}
		for(i=0; i<curr_matrix->get_reaction(j)->num_prod; i++) {
			cnt++;
			row_indices[cnt]=curr_matrix->get_reaction(j)->products[i]->met_num+1;
			col_indices[cnt]=j+1;
			coefficients[cnt]=(double)(curr_matrix->get_reaction(j)->prod_stoich[i]);
		}
		
		
		
	}
	
	
	
}


void Linear_Programmer::set_bounds() 
{
	int i;
	
	for(i=1; i<=curr_matrix->get_num_reactions(); i++) {
		
		switch(curr_matrix->get_reaction(i-1)->type) {
			case BD_FREE:
                //cout<<"Rxn "<<i<<" FREE\n";
				glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
				break;
			case BD_LOWER:
                //cout<<"Rxn "<<i<<" LR: "<<curr_matrix->get_reaction(i-1)->lower_bound<<"\n";
				glp_set_col_bnds(lp, i, GLP_LO, curr_matrix->get_reaction(i-1)->lower_bound, 0.0);
				break;
			case BD_UPPER:
                //cout<<"Rxn "<<i<<" UP: "<<curr_matrix->get_reaction(i-1)->upper_bound<<"\n";
				glp_set_col_bnds(lp, i, GLP_UP, 0.0, curr_matrix->get_reaction(i-1)->upper_bound);
				break;
            case BD_FIXED:
                glp_set_col_bnds(lp, i, GLP_FX, curr_matrix->get_reaction(i-1)->lower_bound, curr_matrix->get_reaction(i-1)->lower_bound);
                break;
			case BD_BOTH:
                //cout<<"Rxn "<<i<<" BOTH: "<<curr_matrix->get_reaction(i-1)->lower_bound<<", "<<curr_matrix->get_reaction(i-1)->upper_bound<<"\n";
				glp_set_col_bnds(lp, i, GLP_DB, curr_matrix->get_reaction(i-1)->lower_bound, curr_matrix->get_reaction(i-1)->upper_bound);
				break;
		}
	}
}


void Linear_Programmer::do_optimization()
{
	do_optimization(false);
}

void Linear_Programmer::do_optimization(bool save_fluxes)
{
	int i, retval;
	double val;
	glp_smcp params;
	glp_init_smcp(&params);
	
	params.presolve=GLP_ON;
    glp_term_out(GLP_OFF);
	glp_load_matrix(lp, non_zero_coeff, row_indices, col_indices, coefficients); 
	glp_simplex(lp, NULL); 
	retval=glp_get_status(lp);
	//cout<<"RETVAL: "<<retval<<endl;
	//success=true;
	if (retval == GLP_OPT) {
		success=true;
		if (save_fluxes == true) {
			for(i=1; i <= curr_matrix->get_num_reactions(); i++) {
				curr_matrix->get_reaction(i-1)->flux =0.0;
				val = glp_get_col_prim(lp, i);
				if (fabs(val) > 0.001)	
						curr_matrix->get_reaction(i-1)->flux = val;
			}
		}
	}
	else success=false;
	
}


double Linear_Programmer::get_opt_val()
{
	return(glp_get_obj_val(lp));
}

void Linear_Programmer::print_soln()
{
	int i;
	double val;
	
	cout<<"Optimized flux through "<<curr_matrix->get_reaction(opt_equation_num)->react_name<<":\t";
	for (i=0; i<curr_matrix->get_reaction(opt_equation_num)->num_react; i++)
		cout<<curr_matrix->get_reaction(opt_equation_num)->react_stoich[i]<<" "<<curr_matrix->get_reaction(opt_equation_num)->reactants[i]->met_name<<" + ";
	cout<<" ---> ";
	
	for (i=0; i<curr_matrix->get_reaction(opt_equation_num)->num_prod; i++)
		cout<<curr_matrix->get_reaction(opt_equation_num)->prod_stoich[i]<<" "<<curr_matrix->get_reaction(opt_equation_num)->products[i]->met_name<<" + ";
	cout<<endl;
	
	val = glp_get_obj_val(lp);
	cout<<"Optimal Solution value: "<<val<<endl;
	
	cout<<"Flux coefficients: \n";
	for(i=1; i <= curr_matrix->get_num_reactions(); i++) {
		curr_matrix->get_reaction(i-1)->flux =0.0;
		val = glp_get_col_prim(lp, i);
		if (fabs(val) > 0.001)	{
            cout<<curr_matrix->get_reaction(i-1)->react_name<<": ";
           /* if (curr_matrix->get_reaction(i-1)->reversible==true)
                cout<<"(R) ";
            else
                cout<<"(I) ";*/
            cout<<val<<endl;
			curr_matrix->get_reaction(i-1)->flux = val;
		}
	}
	
	
	
	//for(i=1; i<6; i++) {
	//	x1=glp_get_row_prim(lp, i);
	//	cout<<"Row variable: "<<i<<": "<<x1<<endl;
	//}
	
	
}


void Linear_Programmer::save_soln()
{
    int i;
    double val;
    
   
    
    val = glp_get_obj_val(lp);
    
    for(i=1; i <= curr_matrix->get_num_reactions(); i++) {
        curr_matrix->get_reaction(i-1)->flux =0.0;
        val = glp_get_col_prim(lp, i);
        if (fabs(val) > 0.001)
            curr_matrix->get_reaction(i-1)->flux = val;
       
    }
    
  
}


void Linear_Programmer::bound_with_current_opt()
{
	int i;
	double val;
	
	
	for(i=1; i <= curr_matrix->get_num_reactions(); i++) {
		curr_matrix->get_reaction(i-1)->flux =0.0;
		val = glp_get_col_prim(lp, i);
		if (fabs(val) > 0.001)	{
						
			if (val >0.0) {
				curr_matrix->get_reaction(i-1)->upper_bound=val;
				curr_matrix->get_reaction(i-1)->type=BD_UPPER;	
			}
			else {
				curr_matrix->get_reaction(i-1)->lower_bound=val;
				curr_matrix->get_reaction(i-1)->type=BD_LOWER;	
			}
		}
	}
	
	set_bounds();
	
}
