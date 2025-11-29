#include "linked_list.cpp"
#include "stoich_mat.h"
#include <string>
#include <algorithm>
#include <cctype>


#define NULL_TOL 1e-12

using namespace std;

//extern std::stringstream dbss;

int Metabolite::operator== (Metabolite test)
{
	int retval=0;
	if(test.met_name ==met_name) retval=1;
	return(retval);
}

Metabolite& Metabolite::operator= (Metabolite &assign_from)
{
	met_name=assign_from.met_name;
	met_num=assign_from.met_num;
	return(*this);
}


Reaction::Reaction()
{
	reactants=products=0;
	react_stoich=prod_stoich=0;
	num_react=num_prod=0;
	reversible=true;
	lethal=false;
	type = BD_FREE;
	lower_bound=upper_bound=flux=0;
}

Reaction& Reaction::operator= (Reaction& assign_from)
{
	int i;
	react_name=assign_from.react_name;
	react_num=assign_from.react_num;
	
	num_react=assign_from.num_react;
	num_prod=assign_from.num_prod;
	
	if (reactants !=0) delete[] reactants;
	if (products != 0) delete[] products;
	if (react_stoich != 0) delete[] react_stoich;
	if (prod_stoich != 0) delete[] prod_stoich;
	
	reversible=assign_from.reversible;
	reactants=0;
	products=0;
	type=assign_from.type;
	upper_bound=assign_from.upper_bound;
	lower_bound=assign_from.lower_bound;
	flux=assign_from.flux;
	exchange=assign_from.exchange;
	
	if (num_react > 0) {
		reactants=new Metabolite*[num_react];
		react_stoich=new double[num_react];
	}
	
	if (num_prod > 0) {
		products=new Metabolite*[num_prod];
		prod_stoich=new double[num_prod];
	}	
	
	for(i=0; i<assign_from.num_react; i++) {
	//	reactants[i]=assign_from.reactants[i];
		react_stoich[i]=assign_from.react_stoich[i];
	}
	
	for(i=0; i<assign_from.num_prod; i++) {
	//	products[i]=assign_from.products[i];
		prod_stoich[i]=assign_from.prod_stoich[i];
	}
	return(*this);
}

void Reaction::knockout_rxn()
{
    old_lower=lower_bound;
    old_upper=upper_bound;
    old_type=type;
    
    upper_bound=0;
    lower_bound=0;
    type = BD_FIXED;
}

void Reaction::restore_rxn()
{
    lower_bound=old_lower;
    upper_bound=old_upper;
    type=old_type;
}

void Reaction::set_products_reactants(std::list<Metabolite*> &reactant_list, std::list<Metabolite*> &product_list, std::list<double> &r_stoich, std::list<double> &p_stoich)
{
	int i;
	
	num_prod=product_list.size();
	num_react=reactant_list.size();
	
	if (reactants !=0) delete[] reactants;
	if (products != 0) delete[] products;
	if (react_stoich != 0) delete[] react_stoich;
	if (prod_stoich != 0) delete[] prod_stoich;
	
	reactants=0;
	products=0;
	
	if (num_react > 0) {
		reactants=new Metabolite*[num_react];
		react_stoich=new double[num_react];
	}
	
	if (num_prod > 0) {
		products=new Metabolite*[num_prod];
		prod_stoich=new double[num_prod];
	}
	

	
    i=0;
    for (std::list<Metabolite*>::const_iterator it = reactant_list.begin(); it != reactant_list.end(); ++it) {
        reactants[i]=(*it);
        i++;
    }
    
    i=0;
    for (std::list<double>::const_iterator it = r_stoich.begin(); it != r_stoich.end(); ++it) {
        react_stoich[i]=(*it);
        i++;
    }
	
    i=0;
    for (std::list<Metabolite*>::const_iterator it = product_list.begin(); it != product_list.end(); ++it) {
        products[i]=(*it);
        i++;
    }
    
    i=0;
    for (std::list<double>::const_iterator it = p_stoich.begin(); it != p_stoich.end(); ++it) {
        prod_stoich[i]=(*it);
        i++;
    }
}


Reaction::~Reaction()
{
	if (reactants !=0) {
		delete[] reactants;
		delete[] react_stoich;
	}
	if(products != 0) {		
		delete[] products;
		delete[] prod_stoich;
	}
}

Reaction_matrix::Reaction_matrix()
{
	num_reactions=num_metabolites=0;
	omit_reaction=-1;
	metabolites=0;
	reactions=0;
	fortran_matrix=0;
	u_mat=0;
	v_mat=0;
	sing_vals=0;
	null_space_span=0;
	
	trans=new char;
#ifdef ___LOCAL_BLAS___
	fortran_cols=new long int;
	fortran_rows=new long int;
	nrhs=new long int;
	lwork=new long int;
	info=new long int;
#else
	fortran_cols=new __CLPK_integer;
	fortran_rows=new __CLPK_integer;
	nrhs=new __CLPK_integer;
	lwork=new __CLPK_integer;
	info=new __CLPK_integer;
#endif
	internal_data=false;
	omit_size=0;
}


Reaction_matrix::Reaction_matrix(string  infile)
{
    //Goal rid this class of incessant iterating through arrays to find metabolite values
    //To be done through a c++ map
    //Metabolite_hash<string metabolite_name, int metabolite_num>
    //At the very end this structure will be iterated once to put all the names
    //and numbers into a single Metabolite class
    
    
   
    fortran_matrix=0;
    sing_vals=0;
    u_mat=0;
    v_mat=0;
    null_space_span=0;
    
    trans=new char;
    
#ifdef ___LOCAL_BLAS___
    fortran_cols=new long int;
    fortran_rows=new long int;
    nrhs=new long int;
    lwork=new long int;
    info=new long int;
#else
    fortran_cols=new __CLPK_integer;
    fortran_rows=new __CLPK_integer;
    nrhs=new __CLPK_integer;
    lwork=new __CLPK_integer;
    info=new __CLPK_integer;
#endif
    internal_data=true;
    omit_reaction=-1;
    omit_size=0;
    
    std::ifstream f_fromfile(infile.c_str());
    
    read_reaction_matrix(f_fromfile);
    f_fromfile.close();
}

Reaction_matrix::Reaction_matrix(std::istream& datass)
{
    //Goal rid this class of incessant iterating through arrays to find metabolite values
    //To be done through a c++ map
    //Metabolite_hash<string metabolite_name, int metabolite_num>
    //At the very end this structure will be iterated once to put all the names
    //and numbers into a single Metabolite class
    
    

    
    fortran_matrix=0;
    sing_vals=0;
    u_mat=0;
    v_mat=0;
    null_space_span=0;
    
    trans=new char;
    
#ifdef ___LOCAL_BLAS___
    fortran_cols=new long int;
    fortran_rows=new long int;
    nrhs=new long int;
    lwork=new long int;
    info=new long int;
#else
    fortran_cols=new __CLPK_integer;
    fortran_rows=new __CLPK_integer;
    nrhs=new __CLPK_integer;
    lwork=new __CLPK_integer;
    info=new __CLPK_integer;
#endif
    internal_data=true;
    omit_reaction=-1;
    omit_size=0;
    
    

    read_reaction_matrix(datass);
    
}



void Reaction_matrix::read_reaction_matrix(std::istream& fin)
{
    int i, j, num_react, num_prod, reaction_cnt, my_hash_pos;
    double stoich;
    char dump;
    string line, new_name;
    bool found;
    //ifstream fin;
    Metabolite a_metabolite, *my_metabolite, *new_metabolite;
    Reaction *a_reaction;
    std:list<double> react_stoich, prod_stoich;
    std::list<Metabolite*> metabolite_list;
    std::list<Metabolite*> reactant_list, product_list;
    std::list<Reaction*> reaction_list;


    react_stoich.clear();
    prod_stoich.clear();
    reactant_list.clear();
    product_list.clear();
    metabolite_list.clear();
    reaction_list.clear();
    
	num_metabolites=0;
	reaction_cnt=0;
	
	Metabolite_hash.clear();
	//cout<<"Initial hash size is "<<Metabolite_hash.size()<<endl;
	getline(fin, line);
	fin>>new_name>>dump;
    //dbss<<"read: "<<new_name;
    
	while(new_name.find("(Exchange") == std::string::npos) {
        
        //cout<<new_name<<"Irr char="<<dump<<endl;
        
        react_stoich.clear();
        prod_stoich.clear();
        reactant_list.clear();
        product_list.clear();
		
		a_reaction=new Reaction;
		
		num_react = num_prod=0;
		a_reaction->react_num=reaction_cnt;
		a_reaction->exchange=false;
		a_reaction->react_name=new_name;
		reaction_cnt++;
		
        
		//cout<<"Reaction: "<<a_reaction->react_name<<endl;
		
		if ((dump == 'i') || (dump == 'I')) a_reaction->reversible=false;
		else a_reaction->reversible=true;
		
		//cout<<"Dump: "<<dump<<endl;
		while (dump != '\n') {
			fin>>stoich>>a_metabolite.met_name;
			//cout<<"Read: |"<<a_metabolite.met_name<<"| s: "<<stoich<<" DUMP: "<<(dump-0)<<endl;
			
			if(Metabolite_hash.find(a_metabolite.met_name) == Metabolite_hash.end()){
				//cout<<"REAC Adding to hash.  Size is now "<<Metabolite_hash.size()<<endl;
                new_metabolite=new Metabolite;
                new_metabolite->met_name=a_metabolite.met_name;
                new_metabolite->met_num=Metabolite_hash.size()-1;
				Metabolite_hash[new_metabolite->met_name] = new_metabolite;
				
				metabolite_list.push_back(new_metabolite);
				
				my_metabolite=new_metabolite;
			}
			else {
				my_metabolite=Metabolite_hash[a_metabolite.met_name];
				//if (strcmp (a_metabolite.met_name, my_metabolite->met_name) != 0)
				//	cout<<"ERROR in matching: "<<a_metabolite.met_name<<" to "<<my_metabolite->met_name<<endl<<flush;
			}

			
			if (stoich <0) {
				reactant_list.push_back(my_metabolite);
				react_stoich.push_back(stoich);
				num_react++;
			}
			else {
				product_list.push_back(my_metabolite);
				prod_stoich.push_back(stoich);
				num_prod++;
			}
			
			
			
			fin.get(dump);
        }
		a_reaction->set_products_reactants(reactant_list, product_list, react_stoich, prod_stoich);
        reaction_list.push_back(a_reaction);
        Reaction_hash[a_reaction->react_name]=a_reaction;
		
		//for(i=0; i<reaction_list.get_list_length(); i++)
		//	cout<<i<<": "<<(*reaction_list.get_nth_element(i)->item())->react_name<<endl;
		
		fin>>new_name>>dump;
        //cout<<"Read "<<new_name<<endl;
	}
	
	getline(fin, line);
    //cout<<"End of line for Exchange "<<line<<endl;
    
    react_stoich.clear();
    prod_stoich.clear();
    reactant_list.clear();
    product_list.clear();
    
	num_react = 0;
	
	while(! fin.eof()) {
		fin>>a_metabolite.met_name>>line;
		
		if (! fin.eof()) {
			
            prod_stoich.clear();
            product_list.clear();
            
			//Get rid of this search and access hash
			if(Metabolite_hash.find(a_metabolite.met_name) == Metabolite_hash.end()){
                new_metabolite=new Metabolite;
                new_metabolite->met_name=a_metabolite.met_name;
                new_metabolite->met_num=Metabolite_hash.size()-1;
                Metabolite_hash[new_metabolite->met_name] = new_metabolite;
                
                metabolite_list.push_back(new_metabolite);
                
                my_metabolite=new_metabolite;
            }
			else{
                my_metabolite=Metabolite_hash[a_metabolite.met_name];
            }
			
						
			product_list.push_back(my_metabolite);
			stoich=1.0;
			prod_stoich.push_back(stoich);
			
			a_reaction=new Reaction;
			
			a_reaction->react_name= "Exchng_" + my_metabolite->met_name;
			a_reaction->exchange=true;
			a_reaction->react_num=reaction_cnt;
			reaction_cnt++;
			
			a_reaction->set_products_reactants(reactant_list, product_list, react_stoich, prod_stoich);
			//reaction_list.add_to_list(a_reaction);
            reaction_list.push_back(a_reaction);
            Reaction_hash[a_reaction->react_name]=a_reaction;
            
			//for(i=0; i<reaction_list.get_list_length(); i++)
			//	cout<<i<<": "<<(*reaction_list.get_nth_element(i)->item())->react_name<<endl;
			
        }
    }

	
	num_reactions=reaction_list.size();
	num_metabolites=metabolite_list.size();
	
	metabolites=new Metabolite*[num_metabolites];
	
    i=0;
    for (std::list<Metabolite*>::const_iterator it = metabolite_list.begin(); it != metabolite_list.end(); ++it) {
        metabolites[i]=(*it);
        metabolites[i]->met_num=i;
        i++;
    }
	
	
	//cout<<"Listing metabolites\n";
	//for(i=0; i<num_metabolites; i++)
	//	cout<<metabolites[i]->met_name<<"\t"<<metabolites[i]->met_num<<endl;
	
	reactions=new Reaction *[num_reactions];
   
    i=0;
    for (std::list<Reaction*>::const_iterator it = reaction_list.begin(); it != reaction_list.end(); ++it) {
        
        reactions[i]=(*it);
        reactions[i]->react_num=i;
        //Reaction_hash[reactions[i]->react_name]=i;
        
        if (reactions[i]->reversible == true)
            reactions[i]->type=BD_FREE;
        
        else
            reactions[i]->type=BD_LOWER;
        
        reactions[i]->lower_bound=0.0;
        reactions[i]->upper_bound=0.0;
        i++;
    }

}


Reaction_matrix& Reaction_matrix::operator= (Reaction_matrix& assign_from)
{
	int i, j, k, target;
	
	if (fortran_matrix != 0) delete[] fortran_matrix;
	if (u_mat != 0) delete[] u_mat;
	if (v_mat != 0) delete[] v_mat;
	if(sing_vals != 0) delete[] sing_vals;
	
	for(i=0; i<num_metabolites; i++) delete metabolites[i];
	for(i=0; i<num_reactions; i++) delete reactions[i];
	
	delete[] metabolites;
	delete[] reactions;
    Metabolite_hash.clear();
    Reaction_hash.clear();
	
	num_metabolites=assign_from.get_num_metabolites();
	num_reactions=assign_from.get_num_reactions();
	
	metabolites=new Metabolite* [num_metabolites];
	reactions=new Reaction* [num_reactions];
	
	
	for(i=0; i<num_metabolites; i++) {
		metabolites[i]=new Metabolite();
		(*metabolites[i])=(*assign_from.get_metabolite(i));
        Metabolite_hash[metabolites[i]->met_name] = metabolites[i];
	}
	
	for(i=0; i<num_reactions; i++) {
		reactions[i]=new Reaction();
		(*reactions[i])=(*assign_from.get_reaction(i));
        Reaction_hash[reactions[i]->react_name]=reactions[i];
        
	}
	
	for(i=0; i<num_reactions; i++) {
		for(j=0; j<reactions[i]->num_react; j++) {
			target=assign_from.get_reaction(i)->reactants[j]->met_num;
			k=0;
			while(metabolites[k]->met_num != target) k++;
			reactions[i]->reactants[j]=metabolites[k];
		}
		for(j=0; j<reactions[i]->num_prod; j++) {
			target=assign_from.get_reaction(i)->products[j]->met_num;
			k=0;
			while(metabolites[k]->met_num != target) k++;
			reactions[i]->products[j]=metabolites[k];
		}
	}
	
	
	return(*this);
}

Metabolite* Reaction_matrix::get_metabolite_by_name(string name)
{
    if (Metabolite_hash.find(name) == Metabolite_hash.end())
        return(0);
    else
        return(Metabolite_hash[name]);
}


Reaction* Reaction_matrix::get_reaction_by_name(string name)
{
    if (Reaction_hash.find(name) == Reaction_hash.end())
        return(0);
    else
        return(Reaction_hash[name]);
}



void Reaction_matrix::set_omit_reaction(int s)
{
	int i, j, new_met_num;
	
	omit_reaction=s;
	omit_size=1;
	
	for(i=0; i<num_metabolites; i++)
		metabolites[i]->visit=false;
	
	new_met_num=0;
	
	for(i=0; i<num_reactions; i++) {
		if (i!= omit_reaction) {
			for(j=0; j<reactions[i]->num_react; j++) {
				if (reactions[i]->reactants[j]->visit == false) {
					reactions[i]->reactants[j]->met_num=new_met_num;
					reactions[i]->reactants[j]->visit=true;
					new_met_num++;
				}
			}
			for(j=0; j<reactions[i]->num_prod; j++) {
				if (reactions[i]->products[j]->visit == false) {
					reactions[i]->products[j]->met_num=new_met_num;
					reactions[i]->products[j]->visit=true;
					new_met_num++;
				}
			}
		}
	}
	
	for(j=0; j<reactions[omit_reaction]->num_react; j++) {
		if (reactions[omit_reaction]->reactants[j]->visit == false) {
			reactions[omit_reaction]->reactants[j]->met_num=new_met_num;
			reactions[omit_reaction]->reactants[j]->visit=true;
			new_met_num++;
		}
	}
	for(j=0; j<reactions[omit_reaction]->num_prod; j++) {
		if (reactions[omit_reaction]->products[j]->visit == false) {
			reactions[omit_reaction]->products[j]->met_num=new_met_num;
			reactions[omit_reaction]->products[j]->visit=true;
			new_met_num++;
		}
	}
}


void Reaction_matrix::real_omit_reaction(int s)
{
	int i,j, cnt, new_met_num;
	Reaction **new_reactions;
	
	omit_reaction=s;
	
	cout<<"Reaction name is "<<reactions[s]->react_name<<endl;
	cout<<reactions[s]->reactants[0]->met_name<<": "<<reactions[s]->reactants[0]->met_num<<endl;
	cout<<reactions[s]->products[0]->met_name<<": "<<reactions[s]->products[0]->met_num<<endl;
	
	new_reactions=new Reaction * [num_reactions-1];
	
	cnt=0;
	for(i=0; i<num_reactions; i++) {
		if (i != omit_reaction) {
			new_reactions[cnt]=reactions[i];
			cnt++;
		}
	}
	
	delete[] reactions;
	reactions=new_reactions;

	omit_reaction=-1;
	num_reactions=num_reactions-1;
	omit_size=0;
	
	new_met_num=0;
	
	for(i=0; i<num_metabolites; i++)
		metabolites[i]->visit=false;
	
//	for(i=0; i<num_metabolites; i++)
//		cout<<metabolites[i]->met_num<<": "<<metabolites[i]->met_num<<" V: "<<metabolite[i]
	
	for(i=0; i<num_reactions; i++) {
		for(j=0; j<reactions[i]->num_react; j++) {
			if (reactions[i]->reactants[j]->visit == false) {
				reactions[i]->reactants[j]->met_num=new_met_num;
				reactions[i]->reactants[j]->visit=true;
				new_met_num++;
			}
		}
		for(j=0; j<reactions[i]->num_prod; j++) {
			if (reactions[i]->products[j]->visit == false) {
				reactions[i]->products[j]->met_num=new_met_num;
				reactions[i]->products[j]->visit=true;
				new_met_num++;
			}
		}
	}
	
	cout<<"Now have "<<new_met_num<<" reactants\n";
}


void Reaction_matrix::create_fortran_stoic_matrix()
{
	int i, j, offset;
	Reaction *dummy;
	
	*fortran_cols=num_reactions-omit_size;
	*fortran_rows=num_metabolites;
	
	//cout<<"Omitsize: "<<omit_size<<" Omit: "<<omit_reaction<<endl;
	//cout<<"M: "<<*fortran_rows<<"\tN:"<<*fortran_cols<<endl;
	
	fortran_matrix=new double [ (*fortran_cols)*(*fortran_rows)];
	
	u_mat=new double [(*fortran_rows)*(*fortran_rows)];
	v_mat=new double [(*fortran_cols)*(*fortran_cols)];
	sing_vals=new double [(*fortran_rows)];
	
	for(i=0; i<*fortran_cols; i++) {
		for(j=0; j<num_metabolites; j++) {
			fortran_matrix[(i*num_metabolites)+j]=0.0;
		}
	}
	
	offset=0;
	for(i=0; i<num_reactions; i++) {
		//cout<<"Reaction is "<<reactions[i]->react_name<<endl;
		if (i != omit_reaction) {
			for(j=0; j<reactions[i]->num_react; j++) {
				//	dummy=reactions[i];
				//	cout<<"Updating "<<reactions[i]->reactants[j]->met_name<<": "<<reactions[i]->reactants[j]->met_num<<endl;
				fortran_matrix[((i-offset)*num_metabolites)+reactions[i]->reactants[j]->met_num]=reactions[i]->react_stoich[j];
			}
			for(j=0; j<reactions[i]->num_prod; j++) {
				fortran_matrix[((i-offset)*num_metabolites)+reactions[i]->products[j]->met_num]=reactions[i]->prod_stoich[j];
			}
		}
		else {offset=1;}
	}
}


void Reaction_matrix::compute_ortho_basis()
{
	int i,j, k, extra_zero_vecs=0, null_cnt, offset;
	double prod;
	bool has_non_null, has_error;
	
	*trans='A';
	
	*lwork=10.0*(*fortran_rows)*(*fortran_rows) +4*(*fortran_rows);
#ifdef ___LOCAL_BLAS___
	iwork = new long int [8*(*fortran_rows)];
#else
	iwork = new __CLPK_integer [8*(*fortran_rows)];
#endif
	work=new double[*lwork];
	
    cout<<"Making fortran_matrix\n";
	if (fortran_matrix ==0)
		create_fortran_stoic_matrix();
	
	//	for(i=0; i<*fortran_cols; i++) {
	//	for(j=0; j<*fortran_rows; j++)
	//		cout<<fortran_matrix[(i*(*fortran_rows))+j]<<"\t";
	//	cout<<"\n";
	//}
	
	
	dgesdd_(trans, fortran_rows, fortran_cols, fortran_matrix, fortran_rows, sing_vals, u_mat, fortran_rows, v_mat, fortran_cols, work, lwork, iwork, info );
	delete[] work;
	
	/*cout<<"Vmatrix\n";
	 for(i=0; i<*fortran_cols; i++) {
	 for(j=0; j<*fortran_cols; j++) 
	 cout<<v_mat[(i*(*fortran_cols))+j]<<"\t";
	 cout<<endl;
	 }
	 
	 cout<<"Umatrix\n";
	 for(i=0; i<*fortran_rows; i++) {
	 for(j=0; j<*fortran_rows; j++) 
	 cout<<u_mat[(i*(*fortran_rows))+j]<<"\t";
	 cout<<endl;
	 }*/
	
	has_error=false;
	
	for(i=0; i<*fortran_rows; i++) {
		prod=0;
		for(j=0; j<(num_reactions-omit_size); j++)
			prod+=v_mat[i+(*fortran_cols)*(j)]*v_mat[i+(*fortran_cols)*(j)];
		if(fabs(prod-1.0) > 1e-7) {
			has_error=true;
			// cerr<<"Error: vector with sv: "<<sing_vals[i]<<" has length: "<<prod<<endl;
		}
	}
	
	if (has_error==true) {
		cout<<"Re-runnning SVD of reaction "<<omit_reaction<<endl;
		if (3*(*fortran_rows)+*fortran_cols > 5*(*fortran_rows)) {*lwork =3*(*fortran_rows)+*fortran_cols;}
		else {*lwork = 5*(*fortran_rows);}
		work = new double [*lwork];
		
		for(i=0; i<*fortran_cols; i++) {
			for(j=0; j<num_metabolites; j++) {
				fortran_matrix[(i*num_metabolites)+j]=0.0;
			}
		}
		
		offset=0;
		for(i=0; i<num_reactions; i++) {
			//cout<<"Reaction is "<<reactions[i]->react_name<<endl;
			if (i != omit_reaction) {
				for(j=0; j<reactions[i]->num_react; j++) {
					//	dummy=reactions[i];
					//	cout<<"Updating "<<reactions[i]->reactants[j]->met_name<<": "<<reactions[i]->reactants[j]->met_num<<endl;
					fortran_matrix[((i-offset)*num_metabolites)+reactions[i]->reactants[j]->met_num]=reactions[i]->react_stoich[j];
				}
				for(j=0; j<reactions[i]->num_prod; j++) {
					fortran_matrix[((i-offset)*num_metabolites)+reactions[i]->products[j]->met_num]=reactions[i]->prod_stoich[j];
				}
			}
			else {offset=1;}
		}		
		dgesvd_(trans, trans, fortran_rows, fortran_cols, fortran_matrix, fortran_rows, sing_vals, u_mat, fortran_rows, v_mat, fortran_cols,
				work, lwork, info );
		delete[] work;
	}
	
	dim_null_space=0;
	has_error=false;
	cout<<"Singular vals:\n";
	for(i=0; i<*fortran_rows; i++) {
		if (fabs(sing_vals[i]) <=NULL_TOL) {
					cout<<sing_vals[i]<<endl;
			dim_null_space++;
			extra_zero_vecs++;
			
			prod=0;
			for(j=0; j<(num_reactions-omit_size); j++)
				prod+=v_mat[i+(*fortran_cols)*(j)]*v_mat[i+(*fortran_cols)*(j)];
			if(fabs(prod-1.0) > 1e-7) {
				has_error=true;
				cerr<<"Error After correction: vector with sv: "<<sing_vals[i]<<" has length: "<<prod<<endl;
			}
		}
		
	}
	
	for(i=num_metabolites; i<(num_reactions-omit_size); i++) {
		has_non_null=false;
		for(j=0; j<(num_reactions-omit_size); j++) {
			if (v_mat[i+(*fortran_cols)*(j)] > NULL_TOL) has_non_null=true;
		}
		
		if (has_non_null == true) dim_null_space++;
	}
	
	cout<<"The null space of this matrix is of dimension: "<<dim_null_space<<endl;
	cout<<"There were "<<extra_zero_vecs<<" extra vectors of singular value 0\n";
	//cout<<"Omit size: "<<omit_size<<endl;
	
	null_space_span=new double*[num_reactions];
	for(i=0; i<num_reactions; i++)
		null_space_span[i]=new double[dim_null_space];
	
	null_cnt=0;
	
	for(i=0; i<num_metabolites; i++) {
		offset=0;
		if  (fabs(sing_vals[i]) <=NULL_TOL) {
			for(j=0; j<num_reactions; j++) {
				if (j != omit_reaction)
					null_space_span[j][null_cnt]=v_mat[i+(*fortran_cols)*(j-offset)]; 
				else {
					offset=1;
					null_space_span[j][null_cnt]=0;
				}
			}
			null_cnt++;
		}
	}
	
	//OFFSET FOR OMIT????	
	
	for(i=num_metabolites; i<(num_reactions-omit_size); i++) {
		
		offset=0;
		has_non_null=false;
		for(j=0; j<(num_reactions-omit_size); j++) {
			if (v_mat[i+(*fortran_cols)*(j)] > NULL_TOL) has_non_null=true;
		}
		if (has_non_null == true) {
			for(j=0; j<num_reactions; j++){
				if (j != omit_reaction)
					null_space_span[j][null_cnt]=v_mat[i+(*fortran_cols)*(j-offset)]; 
				else {
					offset=1;
					null_space_span[j][null_cnt]=0;
				}
			}
			
			null_cnt++;
		}
	}
	
	for (i=0; i<dim_null_space; i++) {
		prod=0;
		for(k=0; k<num_reactions; k++)
			prod+=null_space_span[k][i]*null_space_span[k][i];
		
		if (fabs(1.0-prod) >1e-7) cerr<<"ERROR: vector "<<i<<"not of unit length: "<<prod<<"\n";
		
		for(j=i+1; j<dim_null_space; j++) {
			prod=0;
			for(k=0; k<num_reactions; k++)
				prod+=null_space_span[k][i]*null_space_span[k][j];
			
			if (fabs(prod) >1e-7) cerr<<"ERROR: vectors "<<i<<" and "<<j<<" are not orthogonal: "<<prod<<"\n";
		}
	}
	
	for(i=0; i<num_reactions; i++) {
        cout<<reactions[i]->react_name<<"\t";
        for(j=0; j<dim_null_space; j++) {
            if (fabs(null_space_span[i][j])>1e-10)
                cout<<null_space_span[i][j]<<"\t";
            else
                cout<<"0\t";
        }
		cout<<"\n";
	}
	
	//cout<<"Info: "<<*info<<endl;
	delete[] iwork;
	work=0;
}


Reaction_matrix::~Reaction_matrix()
{
	int i;
	
	if (fortran_matrix !=0) 
		delete[] fortran_matrix;
	
	if (u_mat != 0) delete[] u_mat;
	if (v_mat != 0) delete[] v_mat;
	if (sing_vals != 0) delete[] sing_vals;
	if (null_space_span !=0) {
		for(i=0; i<num_reactions; i++)
			delete[] null_space_span[i];
		delete[] null_space_span;
		
	}
	
	delete fortran_cols;
	delete fortran_rows;
	delete nrhs;
	delete trans;
	delete lwork;
	delete info;
	
	if (internal_data == true) {
		for(i=0; i<num_metabolites; i++)
			delete metabolites[i];
	}
	delete[] metabolites;
	
	if (internal_data == true) {
		for(i=0; i<num_reactions; i++)
			delete reactions[i];
	}
	delete[] reactions;
}



Space_compare::Space_compare()
{
	trans=new char;
	fortran_cols= new long int;
	fortran_rows= new long int;
	nrhs=new long int;
	info=new long int;
	lwork=new long int;
	
	fortran_matrix=fortran_vecs=work=0;
}

double Space_compare::do_comparison (Reaction_matrix *full_matrix, Reaction_matrix *reduced_matrix)
{
	int i, j, start;
	long int k, rows, cols, l;
	char transa ='N', transb='N';
	double *soln_vec, *space_vector, alpha, beta, len, *proj_vec, sqr, retval=0;
	bool all_zero;
	//ofstream fout;
	
	//	fout.open("temp_matrix");
	//	for(i=0; i<reduced_matrix->dim_null_space; i++) {
	//		for(j=0; j<reduced_matrix->get_num_reactions(); j++)
	//			fout<<reduced_matrix->null_space_span[j][i]<<"\t";
	//		fout<<endl;
	//	}
	//	fout.close();
	
	start=0;
	while((start < full_matrix->dim_null_space) && 
		  (fabs(full_matrix->null_space_span[reduced_matrix->get_omit_reaction()][start]) < NULL_TOL)) {start++;}
	
	if (start < full_matrix->dim_null_space) {
		rows=reduced_matrix->dim_null_space;
		cols=1;
		k=reduced_matrix->get_num_reactions();
		fortran_matrix=new double[reduced_matrix->dim_null_space*reduced_matrix->get_num_reactions()];
		
		
		for(i=0; i<reduced_matrix->get_num_reactions(); i++) {
			for(j=0; j<reduced_matrix->dim_null_space; j++) {
				//fortran_matrix[(i*reduced_matrix->dim_null_space)+j]=reduced_matrix->null_space_span[i][j];
				fortran_matrix[(j*reduced_matrix->get_num_reactions())+i]=reduced_matrix->null_space_span[i][j];
				//	cout<<	fortran_matrix[(j*reduced_matrix->get_num_reactions())+i]<<"\t";
			}
			//cout<<endl;
		}
		
		space_vector=new double[full_matrix->get_num_reactions()];
		for(i=0; i<full_matrix->get_num_reactions(); i++) space_vector[i]=full_matrix->null_space_span[i][start];
		//for(i=0; i<full_matrix->get_num_reactions(); i++) cout<<space_vector[i]<<endl;
		soln_vec=new double[reduced_matrix->dim_null_space];
		
		//Multply one of the full null-space vectors by the reduced space
		
		alpha=1.0;
		beta=0.0;
		
		//cblas_dgemm(CblasColMajor,CblasTrans, CblasNoTrans, rows, cols, k, alpha, fortran_matrix, k, space_vector, 
		//			k, beta, soln_vec, rows);
		
#ifdef ___LOCAL_BLAS___
		transa='N';
		transb='N';
		l =1;
		dgemm_ (&transa, &transb, &rows, &cols, &cols, &alpha, fortran_matrix, &k, space_vector, &k, &beta, soln_vec, &rows, &l, &l);	
#else
		cblas_dgemm(CblasColMajor,CblasTrans, CblasNoTrans, rows, cols, k, alpha, fortran_matrix, k, space_vector, 
					k, beta, soln_vec, rows);
#endif	
		//cout<<"Combination to project vector to reduced space\n";
		//for(i=0; i<reduced_matrix->dim_null_space; i++)
		//	cout<<soln_vec[i]<<endl;
		
		proj_vec=new double[reduced_matrix->get_num_reactions()];
	
#if ___LOCAL_BLAS___
		transa='N';
		transb='N';
		rows=reduced_matrix->get_num_reactions();
		cols=1;
		k=reduced_matrix->dim_null_space;
		l=1;
		dgemm_ (&transa, &transb, &rows, &cols, &k, &alpha, fortran_matrix, &rows, soln_vec, &k, &beta, proj_vec, &rows, &l, &l);
#else
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, reduced_matrix->get_num_reactions(), 1, reduced_matrix->dim_null_space,
					alpha, fortran_matrix, reduced_matrix->get_num_reactions(), soln_vec, reduced_matrix->dim_null_space, beta, proj_vec, reduced_matrix->get_num_reactions());
#endif	
		//cout<<"Vector projected to reduced space\n";
		//for(i=0; i<reduced_matrix->get_num_reactions(); i++)
		//	cout<<proj_vec[i]<<endl;
		
		sqr=0;
		for(i=0; i<reduced_matrix->get_num_reactions(); i++) {
			space_vector[i]=full_matrix->null_space_span[i][start]-proj_vec[i];
			sqr+=space_vector[i]*space_vector[i];
		}
		
		//cout<<"Orthologonal vector\n";
		//for(i=0; i<reduced_matrix->get_num_reactions(); i++)
		//	cout<<space_vector[i]<<endl;
		
		sqr=sqrt(sqr);
		//cout<<"Normalizing from length: "<<sqr<<endl;
		for(i=0; i<reduced_matrix->get_num_reactions(); i++) 
			space_vector[i]=space_vector[i]/sqr;
		
		cnt_preclude=0;
		
		for(i=0; i<reduced_matrix->get_num_reactions(); i++) {
			if (space_vector[i] > NULL_TOL) {
				all_zero=true;
				j=0;
				while((j<reduced_matrix->dim_null_space) && (all_zero == true)) {
					if (reduced_matrix->null_space_span[i][j] > NULL_TOL)
						all_zero=false;
					else j++;
				}
				
				if(all_zero == true) cnt_preclude++;
			}
			
			if (space_vector[i] !=0.0) {
				sqr=space_vector[i]*space_vector[i];
				retval+=sqr*(log(sqr)/log(2));
			}
		}
		
		//cout<<"Orthologonal vector\n";
		//for(i=0; i<reduced_matrix->get_num_reactions(); i++)
		//	cout<<space_vector[i]<<endl;
		
		delete[] fortran_matrix;
		fortran_matrix=0;
		delete[] soln_vec;
		delete[] space_vector;
		//dgemm_(transa,transb,fortran_rows,fortran_cols,k,alpha,fortran_matrix,k,space_vector,k,beta,soln_vec,fortran_rows);
		return(retval);
	}
}


Space_compare::~Space_compare()
{
	delete trans;
	delete fortran_cols;
	delete fortran_rows;
	delete nrhs;
	delete info;
	delete lwork;
	
	if (fortran_matrix != 0) delete[] fortran_matrix;
	if(fortran_vecs != 0) delete[] fortran_vecs;
	if(work != 0) delete[] work;
}



RXN_RELATE paired_rxn_line_soln(Reaction_matrix * the_matrix, int rxn1, int rxn2)
{
	int i,start;
	double first_ratio, curr_ratio, len;
	RXN_RELATE ret_val=EMPTY;
	
	i=0;
	
	while ((i<the_matrix->dim_null_space) && ((fabs(the_matrix->null_space_span[rxn1][i]) <= NULL_TOL) || (fabs(the_matrix->null_space_span[rxn2][i]) <= NULL_TOL))) i++;
	
	if (i<the_matrix->dim_null_space) {
		//If we have tranverse the whole null space, then rxn1 and rxn2 are never non-zero together.  Otherwise, we enter this if statement
		ret_val=LINE;
		len=sqrt( (the_matrix->null_space_span[rxn1][i]*the_matrix->null_space_span[rxn1][i]) + 
				 (the_matrix->null_space_span[rxn2][i]*the_matrix->null_space_span[rxn2][i]) );
		first_ratio=fabs(the_matrix->null_space_span[rxn1][i]/len);
		
		start=i;
		for(i=start; i<the_matrix->dim_null_space; i++) {
			if ((fabs(the_matrix->null_space_span[rxn1][i]) > NULL_TOL) && (fabs(the_matrix->null_space_span[rxn2][i]) > NULL_TOL)) {
				len=sqrt( (the_matrix->null_space_span[rxn1][i]*the_matrix->null_space_span[rxn1][i]) + 
						 (the_matrix->null_space_span[rxn2][i]*the_matrix->null_space_span[rxn2][i]) );
				curr_ratio=fabs(the_matrix->null_space_span[rxn1][i]/len);
				
				if (fabs(curr_ratio - first_ratio) > 1e-5) ret_val=PLANE;
			}
			
		}
		
	}
	return(ret_val);
}

void get_bounds (Reaction_matrix *curr_matrix, std::string bound_file)
{
    std::ifstream f_fromfile(bound_file.c_str());
    
    if (f_fromfile.fail()) {
        cerr<<"Error: could not open bound file "<<bound_file<<endl;
        return;
    }
    else {
        get_bounds_read(curr_matrix, f_fromfile);
        f_fromfile.close();
    }
}
void get_bounds (Reaction_matrix *curr_matrix, std::stringstream& bound_ss)
{
    get_bounds_read(curr_matrix, bound_ss);
}
void get_bounds_read (Reaction_matrix *curr_matrix, std::istream& fin)
{
	int i;
	double up, low;
	bool low_inf, high_inf;
	string name, low_string, high_string;
	//ifstream fin;
    Reaction *my_react;
	
	//fin.open(bound_file.c_str());
	
	while(! fin.eof()) {

		fin>>name>>low_string>>high_string;
        std::transform(low_string.begin(), low_string.end(), low_string.begin(), ::tolower);
        std::transform(high_string.begin(), high_string.end(), high_string.begin(), ::tolower);
		//cout<<name<<" Chck by "<<low_string[1]<<"\t"<<high_string[1]<<endl;
        if (low_string.find("inf") == std::string::npos) {
			low_inf=false;
			low=string_to_float(low_string.c_str());
		}
		else {
			low_inf=true;
			low=-1e10;
		}
			
		if (high_string.find("inf") == std::string::npos) {
			high_inf=false;
			up=string_to_float(high_string.c_str());
		}
		else {
			high_inf=true;
			up=1e10;
		}

		//cout<<"SETTING: "<<name<<": "<<low<<"\t"<<up<<endl;
        my_react=curr_matrix->get_reaction_by_name(name);
        
        
		if(my_react == 0) cerr<<"Error: could not find reaction for bound name "<<name<<endl;
		else {
			if (low_inf == false) {
				if (high_inf == false) {
					if (my_react->reversible == true)
						my_react->lower_bound=low;
                    else {
						my_react->lower_bound=fabs(low);
                        if (low > 0.0)
                            cout<<"WARNING: Changing negative bound to positive because "<<my_react->react_name<<" is irreversible\n";
                    }
					my_react->upper_bound=up;
					my_react->type=BD_BOTH;
				}
				else {
					my_react->lower_bound=low;
					my_react->type=BD_LOWER;
				}
			}
			else {
				if (high_inf == false) {
					my_react->upper_bound=up;
					my_react->type=BD_UPPER;
				}
				else {
					my_react->type=BD_FREE;
				}
				
			}
		}
	}
	//fin.close();
}






void get_bounds_var_search (Reaction_matrix *curr_matrix, string bound_file, bool *&bound_var)
{
	int i;
	double up, low;
	bool low_inf, high_inf;
	string name, low_string, high_string, var_string;
	ifstream fin;
    Reaction *my_react;
	
	bound_var=new bool [curr_matrix->get_num_reactions()];
	
	for(i=0; i<curr_matrix->get_num_reactions(); i++) bound_var[i]=false;
	
	fin.open(bound_file.c_str());
	
	if (fin.fail()) {
		cerr<<"Error: could not open bound file "<<bound_file<<endl;
		return;
	}
	
	while(! fin.eof()) {

		fin>>name>>low_string>>high_string>>var_string;
		
        std::transform(low_string.begin(), low_string.end(), low_string.begin(), ::tolower);
        std::transform(high_string.begin(), high_string.end(), high_string.begin(), ::tolower);
        //cout<<name<<" Chck by "<<low_string[1]<<"\t"<<high_string[1]<<endl;
        if (low_string.find("inf") == std::string::npos) {
            low_inf=false;
            low=string_to_float(low_string.c_str());
        }
        else {
            low_inf=true;
            low=-1e10;
        }
        
        if (high_string.find("inf") == std::string::npos) {
            high_inf=false;
            up=string_to_float(high_string.c_str());
        }
        else {
            high_inf=true;
            up=1e10;
        }
	
		//cout<<name<<": "<<low<<"\t"<<up<<endl;
        my_react=curr_matrix->get_reaction_by_name(name);
		
		if(my_react == 0) cerr<<"Error: could not find reaction for bound name "<<name<<endl;
		else {
			
			if (var_string == "VARIABLE") {
				bound_var[my_react->react_num]=true;
			//	cout<<"Allowing bound on "<<curr_matrix->get_reaction(i)->react_name<<" to move\n";
			}
			//else cout<<"Fixing bound on "<<curr_matrix->get_reaction(i)->react_name<<"\n";
			
			if (low_inf == false) {
				if (high_inf == false) {
					if (my_react->reversible == true)
						my_react->lower_bound=low;
					else 
						my_react->lower_bound=fabs(low);
					
					my_react->upper_bound=up;
					my_react->type=BD_BOTH;
				}
				else {
					my_react->lower_bound=low;
					my_react->type=BD_LOWER;
				}
			}
			else {
				if (high_inf == false) {
					my_react->upper_bound=up;
					my_react->type=BD_UPPER;
				}
				else {
					my_react->type=BD_FREE;
				}
				
			}
		}
	}
	fin.close();
}

