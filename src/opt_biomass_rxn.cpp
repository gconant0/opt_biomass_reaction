#include <iostream>
#include <sstream>
#include "lin_program.h"
#include "stoich_mat.h"
#include "gen_dna_funcs.h"



using namespace::std;

#ifdef _FBA_PLOT_
extern int plot_system(double *flux_vals, double *ko_fluxes, string *rxn_names, int num_fluxes, int *flux_cnts, string biomass_name,   string outputfile, std::stringstream *plot_ss, bool streamplot, bool have_ko, bool all_ko);

#endif

int main (int argc, char **argv) {
    int i, j, k, l, opt_num, real_num, new_opt, cnt, non_zero_fluxes=0, flux_cnts[11], loc;
	double upper_bound, old_u[2], old_l[2], val[2], max, inital_bm, *bm_fluxes, *all_fluxes, *flux_vals, *ko_fluxes;
	string react_file, bound_file, limit_rxn, outfile, opt_rxn, knockout_rxn, outplot, *rxn_names;
	ofstream fout;
    stringstream ss;
	bool extra_bound=false, single_kd=false, single_ko=false, double_kd=false,
		lower_bound=false, do_knockout=false, *ran_rxn, same, match, do_plot=false;
	Linear_Programmer *the_lp, *ko_lp;
	Reaction_matrix *the_matrix, *orig_matrix, *new_matrix;
    Reaction *my_react;
	BOUND_TYPE old_type[2];
	
    for(i=0; i<11; i++) flux_cnts[i]=0;
    
	if (argc>=4) {
		react_file=argv[1];
		bound_file=argv[2];
        opt_rxn=argv[3];
        
		opt_num=string_to_int(argv[3]);
		
		if (argc >4) {
            for(i=4; i<argc; i++) {
                if ((argv[i][0] == '-') ||(argv[i][0] == 'Q')) {
                    if ((argv[i][1] == 'k') || (argv[i][1] == 'K')) {
                        knockout_rxn=argv[i];
                        knockout_rxn=knockout_rxn.substr(3, knockout_rxn.length());
                        do_knockout=true;
                    }
                    if (argv[i][1] == 's') {
                        if ((argv[i][3]== 'd') || (argv[i][3]== 'D'))  single_kd=true;
                        else  single_ko=true;
                        
                        outfile=argv[i];
                        outfile=outfile.substr(5, outfile.length);
                    }
                    
                    if ((argv[i][1] == 'd') || (argv[i][1] == 'D')) {
                        double_kd=true;
                        
                        outfile=argv[i];
                        outfile=outfile.substr(5, outfile.length);
                    }
                    if ((argv[i][1] == 'p') || (argv[i][1] == 'P')) {
                        do_plot=true;
                        outplot=argv[i];
                        outplot=outplot.substr(3, outplot.length());
                    }
                }
                else {
                    limit_rxn=argv[i];
                    upper_bound=string_to_float(argv[i+1]);
                    extra_bound=true;
                }
            }
        }
		
		the_matrix=new Reaction_matrix(react_file);
        
     

        my_react=the_matrix->get_reaction_by_name(opt_rxn);
        if (my_react == 0) {
            cerr<<"Error: Cannot find reaction "<<opt_rxn<<" not found in the matrix\n";
            return(-2);
        }
        else opt_num=my_react->react_num;
		
		cout<<"Read "<<the_matrix->get_num_reactions()<<" reactions and "<<the_matrix->get_num_metabolites()<<" metabolites. Optimizing "<<opt_rxn<<" ("<<opt_num<<")\n";;
		
        if (do_knockout==true) {
            my_react=the_matrix->get_reaction_by_name(knockout_rxn);
            //my_react->knockout_rxn();
            if (my_react==0) {
                cerr<<"Error: Cannot find reaction "<<knockout_rxn<<" for knockout\n";
                return(-2);
            }
            else {
                my_react->knockout_rxn();
                cout<<"Knocking out reaction "<<knockout_rxn<<endl;
            }
        }
        
		max =0;
		for(i=0; i<the_matrix->get_reaction(opt_num)->num_react; i++) {
			if(fabs(the_matrix->get_reaction(opt_num)->react_stoich[i]) > max) max=fabs(the_matrix->get_reaction(opt_num)->react_stoich[i]);
		}
		
		if(max > 10000) {
            for(i=0; i<the_matrix->get_reaction(opt_num)->num_react; i++)
                the_matrix->get_reaction(opt_num)->react_stoich[i]=the_matrix->get_reaction(opt_num)->react_stoich[i]/1000;
            
            for(i=0; i<the_matrix->get_reaction(opt_num)->num_prod; i++) 
                the_matrix->get_reaction(opt_num)->prod_stoich[i]=the_matrix->get_reaction(opt_num)->prod_stoich[i]/1000;
        }
		
        
		get_bounds(the_matrix, bound_file);
				
        /*for(i=0; i<the_matrix->get_num_reactions(); i++) {
            cout<<the_matrix->get_reaction(i)->react_name<<": BT "<<the_matrix->get_reaction(i)->type<<" "<< the_matrix->get_reaction(i)->lower_bound<<", "<< the_matrix->get_reaction(i)->upper_bound<<endl;
        }*/
        
        
		the_lp= new Linear_Programmer(the_matrix, opt_num);
		the_lp->do_optimization();
		the_lp->print_soln();
        inital_bm=the_lp->get_opt_val();
        
		for(i=0; i<the_matrix->get_num_reactions(); i++) {
            the_matrix->get_reaction(i)->base_flux=the_matrix->get_reaction(i)->flux;
			the_matrix->get_reaction(i)->min=(the_matrix->get_reaction(i)->flux/the_lp->get_opt_val());
			the_matrix->get_reaction(i)->max=(the_matrix->get_reaction(i)->flux/the_lp->get_opt_val());
		}
		
#ifdef _FBA_PLOT_
            
        if (do_knockout==false) {
                for (i=0; i<the_matrix->get_num_reactions(); i++) {
                    if (the_matrix->get_reaction(i)->flux != 0.0) non_zero_fluxes++;
                }
                flux_vals=new double [non_zero_fluxes];
                rxn_names=new string [non_zero_fluxes];
                
                cnt=0;
        
            
                
                for (i=0; i<the_matrix->get_num_reactions(); i++) {
                    if (the_matrix->get_reaction(i)->flux != 0.0) {
                        flux_vals[cnt]=the_matrix->get_reaction(i)->flux;
                        rxn_names[cnt]=the_matrix->get_reaction(i)->react_name;
                        cnt++;
                    }
                }
            }
            else {
                all_fluxes=new double [the_matrix->get_num_reactions()];
                for (i=0; i<the_matrix->get_num_reactions(); i++) {
                    all_fluxes[i]=the_matrix->get_reaction(i)->flux;
                }
                my_react->restore_rxn();
                the_lp->set_bounds();
                the_lp->do_optimization();
                the_lp->save_soln();
                the_lp->print_soln();
                non_zero_fluxes=0;
                for (i=0; i<the_matrix->get_num_reactions(); i++) {
                    if ((the_matrix->get_reaction(i)->flux != 0.0) ||(all_fluxes[i] !=0)) non_zero_fluxes++;
                }
                flux_vals=new double [non_zero_fluxes];
                rxn_names=new string [non_zero_fluxes];
                ko_fluxes=new double [non_zero_fluxes];
                
                cnt=0;
                for (i=0; i<the_matrix->get_num_reactions(); i++) {
                    if ((the_matrix->get_reaction(i)->flux != 0.0) ||(all_fluxes[i] !=0)) {
                        flux_vals[cnt]=the_matrix->get_reaction(i)->flux;
                        ko_fluxes[cnt]=all_fluxes[i];
                        rxn_names[cnt]=the_matrix->get_reaction(i)->react_name;
                        cnt++;
                    }
                }
                
            }
            
            
#endif
        
		
		if (extra_bound == true) {
            my_react=the_matrix->get_reaction_by_name(limit_rxn);
            
			if(my_react == 0) cerr<<"Error: could not find reaction for bound name "<<limit_rxn<<endl;
			else {
				if (upper_bound > 0) {
					my_react->upper_bound=upper_bound;
					my_react->type=BD_UPPER;
				}
				else {
					my_react->lower_bound=upper_bound;
					my_react->type=BD_LOWER;
				}
				
			}
			
			the_lp->set_bounds();
			the_lp->do_optimization();
			the_lp->print_soln();
            

			
		}
		
		if (single_kd ==true) {
			orig_matrix=new Reaction_matrix();
			(*orig_matrix)=(*the_matrix);
			cout<<"Beginning single knockdowns\n";
			fout.open(outfile.c_str());
			for (i=0; i<the_matrix->get_num_reactions(); i++) {				
				if (the_matrix->get_reaction(i)->flux != 0.0) {
					if (the_matrix->get_reaction(i)->flux > 0) {
						the_matrix->get_reaction(i)->upper_bound=the_matrix->get_reaction(i)->flux/2.0;
						the_matrix->get_reaction(i)->type=BD_UPPER;
					}
					else {
						the_matrix->get_reaction(i)->lower_bound=the_matrix->get_reaction(i)->flux/2.0;
						the_matrix->get_reaction(i)->type=BD_LOWER;
					}
					the_lp->set_bounds();
					the_lp->do_optimization();
					if (the_lp->run_success() == true) {
						val[0]=the_matrix->get_reaction(i)->flux/2.0;

						fout<<the_matrix->get_reaction(i)->react_name<<"\t"<<val[0]<<"\t"<<the_lp->get_opt_val()<<endl;
						cout<<the_matrix->get_reaction(i)->react_name<<"\t"<<val[0]<<"\t"<<the_lp->get_opt_val()<<endl;
                    
					}
					
					(*the_matrix)=(*orig_matrix);
					the_lp->set_bounds();
					the_lp->do_optimization();
				}
			}
			fout.close();
		}
		if (single_ko == true) {
			ran_rxn=new bool [the_matrix->get_num_reactions()];
            bm_fluxes=new double [the_matrix->get_num_reactions()];
            for (i=0; i<the_matrix->get_num_reactions(); i++) {
                ran_rxn[i]=false;
                bm_fluxes[i]=0.0;
            }
			cout<<"Beginning single knockouts\n";
			fout.open(outfile.c_str());
			for (i=0; i<the_matrix->get_num_reactions(); i++) {				
				if ((i != opt_num) && (the_matrix->get_reaction(i)->exchange == false)) {
					if (the_matrix->get_reaction(i)->flux != 0.0) {

                        the_matrix->get_reaction(i)->knockout_rxn();
                        the_lp->set_bounds();
                        the_lp->do_optimization();
                        if ((the_lp->run_success() == true) && (the_lp->get_opt_val() >1e-4)) {
                            ran_rxn[i]=true;
                            bm_fluxes[i]=the_lp->get_opt_val();
                            the_lp->print_soln();
                            
                            loc=1;
                            while((loc<11)&&
                                  (((double)loc/10)<= (bm_fluxes[i]/inital_bm))) loc++;
                            flux_cnts[loc]++;
                            
                            for(j=0; j<the_matrix->get_num_reactions(); j++) {
                                
                                
                                if (the_matrix->get_reaction(j)->min > (the_matrix->get_reaction(j)->flux/the_lp->get_opt_val()))
                                    the_matrix->get_reaction(j)->min=(the_matrix->get_reaction(j)->flux/the_lp->get_opt_val());
                                if (the_matrix->get_reaction(j)->max < (the_matrix->get_reaction(j)->flux/the_lp->get_opt_val()))
                                    the_matrix->get_reaction(j)->max=(the_matrix->get_reaction(j)->flux/the_lp->get_opt_val());
                            }
                        }
                        else {
                            ran_rxn[i]=true;
                            flux_cnts[0]++;
                            the_matrix->get_reaction(i)->lethal=true;
                        }
                        the_matrix->get_reaction(i)->restore_rxn();
                    }
                    
                }
            }
			
            fout<<"RxnID\tKnockedOut\tLethal\tRatioBMProd\tMin\tMax\tMaxMag\n";
			for (i=0; i<the_matrix->get_num_reactions(); i++) {
                if (the_matrix->get_reaction(i)->flux != 0.0) 
					cout<<the_matrix->get_reaction(i)->react_name<<"\t"<<the_matrix->get_reaction(i)->lethal<<"\n";
                fout<<the_matrix->get_reaction(i)->react_name<<"\t";
                if (ran_rxn[i]==true) fout<<"Yes\t";
                else fout<<"No\t";
                if (the_matrix->get_reaction(i)->lethal==true) fout<<"Yes\t";
                else fout<<"No\t";
                
                if (ran_rxn[i]==true)
                    fout<<bm_fluxes[i]/inital_bm<<"\t";
                else
                    fout<<"NA\t";
                fout<<the_matrix->get_reaction(i)->min<<"\t"<<the_matrix->get_reaction(i)->max<<"\t";
				if (fabs(the_matrix->get_reaction(i)->min) > fabs(the_matrix->get_reaction(i)->max)) fout<<fabs(the_matrix->get_reaction(i)->min)<<endl;
				else fout<<fabs(the_matrix->get_reaction(i)->max)<<endl;
			}
			
			
			delete[] ran_rxn;
			fout.close();
		}
		if (double_kd== true) {
			fout.open(outfile.c_str());
			for (i=0; i<the_matrix->get_num_reactions(); i++) {
				for(j=i+1; j<the_matrix->get_num_reactions(); j++) {
					if ((the_matrix->get_reaction(i)->flux != 0.0) && (the_matrix->get_reaction(j)->flux != 0.0)) {
						old_u[0]=the_matrix->get_reaction(i)->upper_bound;
						old_l[0]=the_matrix->get_reaction(i)->lower_bound;
						old_type[0]=the_matrix->get_reaction(i)->type;
						
						old_u[1]=the_matrix->get_reaction(j)->upper_bound;
						old_l[1]=the_matrix->get_reaction(j)->lower_bound;
						old_type[1]=the_matrix->get_reaction(j)->type;
						
						if (the_matrix->get_reaction(i)->flux > 0) {
							the_matrix->get_reaction(i)->upper_bound=the_matrix->get_reaction(i)->flux/2.0;
							the_matrix->get_reaction(i)->type=BD_UPPER;
						}
						else {
							the_matrix->get_reaction(i)->lower_bound=the_matrix->get_reaction(i)->flux/2.0;
							the_matrix->get_reaction(i)->type=BD_LOWER;
						}
						
						if (the_matrix->get_reaction(j)->flux > 0) {
							the_matrix->get_reaction(j)->upper_bound=the_matrix->get_reaction(j)->flux/2.0;
							the_matrix->get_reaction(j)->type=BD_UPPER;
						}
						else {
							the_matrix->get_reaction(j)->lower_bound=the_matrix->get_reaction(j)->flux/2.0;
							the_matrix->get_reaction(j)->type=BD_LOWER;
						}
						
						the_lp->set_bounds();
						the_lp->do_optimization();
						if (the_lp->run_success() == true) {

							val[0]=the_matrix->get_reaction(i)->flux/2.0;
							val[1]=the_matrix->get_reaction(j)->flux/2.0;
							fout<<the_matrix->get_reaction(i)->react_name<<"\t"<<val[0]<<"\t"<<the_matrix->get_reaction(j)->react_name
								<<"\t"<<val[1]<<"\t"<<the_lp->get_opt_val()<<endl;
						
							cout<<i<<"\t"<<j<<the_matrix->get_reaction(i)->react_name<<"\t"<<val[0]<<"\t"<<the_matrix->get_reaction(j)->react_name
								<<"\t"<<val[1]<<"\t"<<the_lp->get_opt_val()<<endl;
						}
						the_matrix->get_reaction(i)->upper_bound=old_u[0];
						the_matrix->get_reaction(i)->lower_bound=old_l[0];
						the_matrix->get_reaction(i)->type=old_type[0];
						the_matrix->get_reaction(j)->upper_bound=old_u[1];
						the_matrix->get_reaction(j)->lower_bound=old_l[1];
						the_matrix->get_reaction(j)->type=old_type[1];
						the_lp->set_bounds();
						the_lp->do_optimization();
					}
				}
			}
			
			fout.close();
		}
		
#ifdef _FBA_PLOT_
        if (do_plot == true) {
            if (single_ko==true) do_knockout=false;
            cout<<"Plotting: "<<do_knockout<<" and "<<single_ko<<" with "<<non_zero_fluxes<<" reactions with flux"<<endl;
            plot_system(flux_vals, ko_fluxes, rxn_names, non_zero_fluxes, flux_cnts, opt_rxn,   "test.ps", &ss, false,  do_knockout, single_ko);
        }
#endif
		
		return(0);
	}
	else {
        cerr<<"Usage: opt_biomass_rxn <matrix file> <boundary file> Biomass_Rxn_Name ((constrain_Rxn) (constraint val)) (-k:RxnName) (-skd:<output file>)  (-sko:<output file>) (-dkd:<output file>) (-p:<plotfile>)\n";
	}

}
