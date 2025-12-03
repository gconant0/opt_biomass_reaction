#include <cstdlib>
#include <cstdio>
#include <plot.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <cstring>


using namespace::std;

#define NUM_BLOCKS 4

extern int plot_system(double *flux_vals, double *ko_fluxes, string *rxn_names, int num_fluxes, int *flux_cnts, string biomass_name,   string outputfile, std::stringstream *plot_ss, bool streamplot, bool have_ko, bool all_ko)
{
    int realx, my_steps=1, frame_steps, freeze, s, realy, i, j, k,  thandle, intlog, intlow, digit, ymarks=10, xmarks=10, maxcnt, my_block;
    double spacex, spacey, border, xrange, yrange, max_flux, xcenter, ycenter, font_size,x1,x2,y1,y2, max_x, spacer, zero_xoff,
        y_topl, x_topl, logval, newval, ymark_l, mark, xmark_l, zero_off, box_size, box_area, val, val2, offset, min_flux,
        line_off, dot_size, max_label, plot_area, block_starts[NUM_BLOCKS], label_num, old_font;
    string *my_names;
    char* buf_pt;
    char  attrib[100], value[100], size_string[25], asteps[50], delay[10];
    FILE *outfile;
    string colors[6], label, print_name;
    std::size_t mypos, size, found;
    
    colors[0]="midnightblue";
    colors[4]="red";
    colors[1]="forestgreen";
    colors[2]="mediumslateblue";
    colors[3]="greenyellow";
    colors[5]="lightslategrey";
    
    max_flux=fabs(flux_vals[0]);
    min_flux=-1;
    
    my_names=new string[num_fluxes];
    
    for(i=0; i<num_fluxes; i++) {
        if (fabs(flux_vals[i])>max_flux) max_flux=fabs(flux_vals[i]);
        if (flux_vals[i] != 0) {
            if ((min_flux <0) ||(fabs(flux_vals[i])<min_flux)) min_flux=fabs(flux_vals[i]);
        }
        if (have_ko == true) {
            if (fabs(ko_fluxes[i])>max_flux) max_flux=fabs(ko_fluxes[i]);
            if (ko_fluxes[i] != 0) {
                if ((min_flux <0) ||(fabs(ko_fluxes[i])<min_flux)) min_flux=fabs(ko_fluxes[i]);
            }
        }
        print_name=rxn_names[i];
        if(print_name.find("Exchng") != std::string::npos) {
            print_name=print_name.substr(found+6, print_name.length());
            print_name= "Ex" + print_name;
            
            if(print_name.find("_e") != std::string::npos) {
                print_name=print_name.substr(0, print_name.length()-2);
            }
            
            
        }
        
        if(print_name.find("_biomass") != std::string::npos) {
           // cout<<"Found "<<print_name<<" at "<<found;
            print_name=print_name.substr(0, print_name.length()-8);
           // cout<<" to "<<print_name<<endl;
        }
        
        if(print_name.find("_ext_b") != std::string::npos) {
           // cout<<"Found "<<print_name<<" at "<<found;
            print_name=print_name.substr(0, print_name.length()-6);
           // cout<<" to "<<print_name<<endl;
        }
        if(print_name.find("_cyt_b") != std::string::npos) {
            print_name=print_name.substr(0, print_name.length()-6);
        }
        
        
        my_names[i]=print_name;
        
        
    }
    
   
    //cout<<"Plotting "<<num_fluxes<<" flux values (max is "<<max_flux<<" and min is "<<min_flux<<") to test.ps\n";
    max_flux=max_flux*1.05;

    logval=log(max_flux)/log(10.0);
    intlog=(int)logval+1;
    
    min_flux = min_flux*0.9;
    logval=log(min_flux)/log(10.0);
    intlow=(int)logval;
    if ((double)intlow > logval) intlow--;
    //xmark_l=intlog;
    xmarks=intlog-intlow;
    //cout<<"Int log is "<<intlog<<endl;
    //xmark_l=newval/(double)xmarks;
  
   //cout<<"Log val is "<<logval<<" and int is "<<intlog<<" and lgval low is "<<logval<<" intlow is "<<intlow<<endl;
    
    
    spacex=500.0;
    spacey=500.0;
    
    if(have_ko==true)
        line_off=15;
    else
        line_off=0;
    
    if ((max_flux>10000) ||(min_flux <1))  {
        border=105;
        spacer=20;
    }
    else {
        border=85;
        spacer=12;
    }
    zero_off=(spacey-(border+spacer))*0.015;
    
    xrange=(spacex-(border+spacer+line_off));
    yrange=(spacey-(border+spacer+line_off));
    
    box_area = yrange/(((double)(num_fluxes+NUM_BLOCKS/2)/(double)NUM_BLOCKS));
    box_size=box_area*0.8;
    
    //cout<<"Yrange is "<<yrange<<" NF: "<<num_fluxes<<" with bs and ba: "<<box_size<<" ... "<<box_area<<endl;
    
    if (streamplot == false)
        outfile=fopen(outputfile.c_str(), "w");
    else {
        size=0;
        //printf ("At open buf = `%s', size = %zu\n", buf_pt, size);
        outfile= open_memstream (&buf_pt, &size);
        
        //fprintf(outfile, "TESTING");
        fflush(outfile);
        //printf ("At test buf = `%s', size = %zu\n", buf_pt, size);
    }
    
    
    if (streamplot == false) {
        freeze=0;
        strcpy(attrib, "PAGESIZE");
        strcpy(value, "letter,xsize=5in,ysize=5in");
        pl_parampl (attrib, value);
        /* create a Postscript Plotter that writes to standard output */
        
        if ((thandle = pl_newpl ("ps", stdin, outfile, stderr)) < 0) {
              fprintf (stderr, "Couldn't create Plotter\n");
              return 1;
        }
       
        pl_selectpl (thandle);       /* select the Plotter for use */
        dot_size=2.1;
    }
    else {
            strcpy(size_string, "1900x1900");
            realx=1900;
            realy=1900;
      
 
        /* set a Plotter parameter */
        pl_parampl ("BITMAPSIZE", size_string);
        
        if ((thandle = pl_newpl ("gif", stdin, outfile, stderr)) < 0) {
            fprintf (stderr, "Couldn't create Plotter\n");
            return 1;
        }
       
        pl_selectpl (thandle);       /* select the Plotter for use */
        dot_size=1.7;
    }

    
    
    if (pl_openpl () < 0)       /* open Plotter */ {
      fprintf (stderr, "Couldn't open Plotter\n");
      return 1;
    }
    
    if ((have_ko==false) &&(all_ko==false)){
        if (streamplot == true) {
            pl_fontname ("HersheySans-BoldOblique"); /* choose a Postscript font */
            pl_ftextangle (0);         /* text inclination angle (degrees) */
            
        }
        else {
            pl_fontname("Helvetica-Oblique");
            pl_ftextangle (0);         /* text inclination angle (degrees) */
            
        }
        if (num_fluxes>400)
            pl_ffontsize (4.4);
        else
            pl_ffontsize (6);
        
        max_label=pl_flabelwidth (my_names[0].c_str());
        for(i=1; i<num_fluxes; i++) {
        
            if (pl_flabelwidth (my_names[i].c_str())>max_label)
                max_label=pl_flabelwidth (my_names[i].c_str());
        }
        max_label=max_label*1.08;
        plot_area = ((spacex-(border+line_off+spacer))-(((double)NUM_BLOCKS-1)*max_label))/((double)NUM_BLOCKS);
        //cout<<"Max le "<<max_label<<" PA: "<<plot_area<<" SP: "<<spacex<<endl;
    }
    
    
    pl_fspace (0.0, 0.0, spacex, spacey); /* specify user coor system */
    pl_flinewidth (1.0);       /* line thickness in user coordinates */
    pl_pencolorname ("black");    /* path will be drawn in red */
    
   
        
        
    if (streamplot == true) {
        pl_fontname ("HersheySans-BoldOblique"); /* choose a Postscript font */
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        pl_ffontsize (14);
        font_size=17.0;
    }
    else {
        pl_fontname("Helvetica-Oblique");
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        pl_ffontsize (14);
        font_size=17;
    }
    
    pl_flinewidth (0.5);
    if (have_ko==false) {
        for(k=0; k<NUM_BLOCKS; k++) {
            
            x1=border+line_off+((double)k*(max_label+plot_area));
            
            block_starts[k]=x1;
            //cout<<"Block "<<k<<" x start: "<<block_starts[k]<<endl;
            x2=x1+plot_area;
            pl_fline(x1, border+line_off, x2, border+line_off);
            pl_fline(x1, spacey-spacer, x2, spacey-spacer);
            pl_fline(x2, border+line_off, x2, spacey-spacer);
            pl_fline(x1, border+line_off, x1, spacey-spacer);
        }
        
       
    }
    else {
        pl_fline(border+line_off, border+line_off, spacex-spacer, border+line_off);
        pl_fline(border+line_off, spacey-spacer, spacex-spacer, spacey-spacer);
        pl_fline(spacex-spacer, border+line_off, spacex-spacer, spacey-spacer);
        pl_fline(border+line_off, border+line_off, border+line_off, spacey-spacer);
       
    }
    
    if (all_ko==false) {
        if(have_ko==true) {
            pl_fline(border+line_off, border, spacex-spacer, border);
            pl_fline(border, border+line_off, border, spacey-spacer);
        }
        
        
        if (have_ko==false) {
            pl_fmove(0.4*spacex+border, 0.6*border);
            pl_alabel('c', 'c', "Flux (arbitrary units)");
        }
        else {
            pl_fmove(0.4*spacex+border, 0.6*border);
            pl_alabel('c', 'c', "Wild-type flux (arbitrary units)");
            
            pl_fmove(0.5*border, 0.4*spacey+border);
            pl_ftextangle (90);
            pl_alabel('c', 'c', "Knockout flux (arbitrary units)");
            
        }
    }
    else {
        pl_fmove(0.4*spacex+border, 0.6*border);
        pl_alabel('c', 'c', "Relative flux (KO/WildType)");
        
        pl_fmove(0.5*border, 0.4*spacey+border);
        pl_ftextangle (90);
        pl_alabel('c', 'c', "Number of Reactions");
    }
    
    pl_ftextangle(0.0);

    if (streamplot == true) {
        pl_fontname ("HersheySans-BoldOblique"); /* choose a Postscript font */
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        if (max_flux>10000)
            old_font=10;
            
        else
            old_font=12;
        font_size=17.0;
    }
    else {
        pl_fontname("Helvetica-Oblique");
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        old_font=10;
        font_size=17;
    }
    
    pl_ffontsize (old_font);
    
    if (all_ko==false) {
    
        if (have_ko==true) {
            for(i=0; i<=xmarks; i++) {
                
                // cout<<"Printing "<<label<<" at "<<newval<<" i "<<i<<" m "<<ymark_l<<" mark v "<<mark<<endl;
                newval=((double)(i))*(xrange/(double)(intlog-intlow))+(border+line_off);
                
                pl_fmove(newval, 0.88*border);
                mark=(double)i;
                std::stringstream label_ss2;
                
                label_num=pow(10, i+intlow);
                
                if ((label_num<100000) && (label_num>0.0001)) {
                    label_ss2<<label_num;
                    pl_alabel('c', 'c', label_ss2.str().c_str());
                }
                else {
                    pl_alabel('c', 'c', "10");
                    
                    pl_ffontsize (old_font-2);
                    pl_fmove(newval+1.05*(pl_flabelwidth ("10")), 0.9*border);
                    //if (label_num<1)label_ss2<<"-";
                    label_ss2<<i+intlow;
                    pl_alabel('c', 'c', label_ss2.str().c_str());
                    
                    pl_ffontsize (old_font);
                    
                }
                
                //label_ss2<<pow(10, i+intlow);
                //pl_alabel('c', 'c', label_ss2.str().c_str());
                pl_fline(newval, border+line_off, newval, 1.05*(border+line_off));
            }
            for(i=0; i<=xmarks; i++) {
                // cout<<"Printing "<<label<<" at "<<newval<<" i "<<i<<" m "<<ymark_l<<" mark v "<<mark<<endl;
                newval=((double)(i))*(xrange/(double)(intlog-intlow))+(border+line_off);
                
                
                mark=(double)i;
                std::stringstream label_ss2;
                
                label_num=pow(10, i+intlow);
                
                if ((label_num<100000) && (label_num>0.0001)) {
                    label_ss2<<label_num;
                    pl_fmove(0.88*border-0.35*(pl_flabelwidth (label_ss2.str().c_str())), newval);
                    pl_alabel('c', 'c', label_ss2.str().c_str());
                    //pl_alabel('c', 'c', label_ss2.str().c_str());
                }
                else {
                    pl_fmove(0.88*border-0.6*(pl_flabelwidth ("10")), newval);
                    pl_alabel('c', 'c', "10");
                    pl_ffontsize (old_font-2);
                    if (label_num<1)
                        pl_fmove(0.88*border+0.4*(pl_flabelwidth ("10")), newval+0.01*yrange);
                    else
                        pl_fmove(0.88*border+0.3*(pl_flabelwidth ("10")), newval+0.01*yrange);
                
                    label_ss2<<i+intlow;
                    pl_alabel('c', 'c', label_ss2.str().c_str());
                    
                    pl_ffontsize (old_font);
                    
                }
                
                label_ss2<<pow(10,i+intlow);
                //pl_ffontsize (8);
                //pl_fmove(0.88*border-0.35*(pl_flabelwidth (label_ss2.str().c_str())), newval);
                //pl_alabel('c', 'c', label_ss2.str().c_str());
                pl_fline(border+line_off, newval, 1.05*(border+line_off), newval);
            }
        }
        else {
           
            
            pl_ffontsize (old_font-1);
            
            for(k=0; k<NUM_BLOCKS; k++) {
                for(i=0; i<=xmarks; i++) {
                    
                    // cout<<"Printing "<<label<<" at "<<newval<<" i "<<i<<" m "<<ymark_l<<" mark v "<<mark<<endl;
                    newval=((double)(i))*(plot_area/(double)(intlog-intlow))+block_starts[k];
                    
                    pl_fmove(newval, 0.88*border);
                    mark=(double)i;
                    std::stringstream label_ss2;
                    
                    if ((i%3)==0) {
                        label_num=pow(10, i+intlow);
                        
                        if (label_num<100000) {
                            label_ss2<<label_num;
                            pl_alabel('c', 'c', label_ss2.str().c_str());
                        }
                        else {
                            pl_alabel('c', 'c', "10");
                            
                            pl_ffontsize (old_font-2);
                            pl_fmove(newval+pl_flabelwidth ("10"), 0.9*border);
                            label_ss2<<i+intlow;
                            pl_alabel('c', 'c', label_ss2.str().c_str());
                            
                            pl_ffontsize (old_font);
                            
                        }
                    }
                    pl_fline(newval, border+line_off, newval, 1.05*(border+line_off));
                    
                }
            }
            pl_ffontsize (old_font);
            
        }
        
        if (streamplot == true) {
            old_font=9;
  
        }
       
        else
            old_font=7;
        
        pl_filltype(0);
        
        pl_flinewidth (1.0);
        
        pl_pencolorname (colors[0].c_str());
        
        
        
        y1=border;
        
        pl_flinewidth (0.0);
        pl_filltype(1);
        
        if (have_ko==false) {
            pl_filltype(0);
            pl_pencolorname ("black");
            pl_flinewidth (0.5);
            
            pl_ffontsize (8);
            
            pl_fbox (0.2*border, 0.6*border, 0.8*border, 0.9*border);
            pl_filltype(1);
            for(i=0; i<2; i++) {
                pl_fillcolorname(colors[i].c_str());
                pl_fbox (0.25*border, i*10.0+0.65*border, 0.45*border,8+i*10.0+0.65*border);
                
                pl_fmove(0.6*border, i*10.0+0.67*border);
                
                if (i==0) {
                    pl_alabel('c', 'c', ">0");
                }
                if (i==1) {
                    pl_alabel('c', 'c', "<0");
                }
            }
            
            if (num_fluxes>400)
                pl_ffontsize (4.4);
            else
                pl_ffontsize (6);
            
            offset=pl_flabelwidth (my_names[0].c_str());
            for(i=1; i<num_fluxes; i++) {
                if (pl_flabelwidth (my_names[i].c_str())>offset)
                    offset=pl_flabelwidth (my_names[i].c_str());
            }
            
            
            for(i=0; i<num_fluxes; i++) {
                //pl_pencolorname ("white");
                if (rxn_names[i]==biomass_name) {
                    pl_fillcolorname(colors[4].c_str());
                    
                }
                else {
                    if (flux_vals[i]>0)
                        pl_fillcolorname(colors[0].c_str());
                    else
                        pl_fillcolorname(colors[1].c_str());
                }
                my_block=i%NUM_BLOCKS;
                
                val =((log(fabs(flux_vals[i]))/log(10))-intlow)/(intlog-intlow);
                //cout<<"Log of flux for "<<rxn_names[i]<<" is "<<log(fabs(flux_vals[i]))/log(10)<<endl;
                //pl_fbox (border, y1, (double)((fabs(flux_vals[i])/max_flux)*xrange)+border, y1+box_size);
                pl_fbox (block_starts[my_block], y1, val*plot_area+block_starts[my_block], y1+box_size);
                
                
                pl_fmove(block_starts[my_block]-(0.59*pl_flabelwidth (my_names[i].c_str())), y1+(0.4*box_size));
                pl_ftextangle(0.0);
                pl_alabel ('c', 'c',my_names[i].c_str());
                if (my_block==(NUM_BLOCKS-1))
                    y1+=box_area;
                
            }
        }
        else {
           
            
            
            pl_filltype(0);
            pl_pencolorname ("black");
            pl_flinewidth (0.5);
            
            
            pl_ffontsize (8);
            
            pl_fbox (0.36*border, 0.55*border, 0.96*border, 0.96*border);
            pl_filltype(1);
            for(i=0; i<5; i++) {
                pl_fillcolorname(colors[i].c_str());
                pl_fcircle (0.42*border, i*8.0+0.6*border, 2);
                
                pl_fmove(0.69*border, i*8.0+0.59*border);
                
                if (i==0) {
                    pl_alabel('c', 'c', "Both >0");
                }
                if (i==2) {
                    pl_alabel('c', 'c', "wt>0, ko<0");
                }
                if (i==1) {
                    pl_alabel('c', 'c', "Both <0");
                }
                if (i==3) {
                    pl_alabel('c', 'c', "wt<0, ko>0");
                }
                if (i==4) {
                    pl_alabel('c', 'c', "Biomass");
                }
            }
            
            
            pl_flinewidth (0.25);
            pl_pencolorname (colors[5].c_str());
            pl_fline(border+line_off, border+line_off, spacex-spacer, spacey-spacer);
            pl_flinewidth (0.0);
            pl_pencolorname ("black");
            
            for(i=0; i<num_fluxes; i++) {
                //pl_pencolorname ("white");
                if (rxn_names[i]==biomass_name) {
                    pl_fillcolorname(colors[4].c_str());
                    
                }
                else {
                    if (flux_vals[i]>0) {
                        if (ko_fluxes[i] >0)
                            pl_fillcolorname(colors[0].c_str());
                        else
                            pl_fillcolorname(colors[2].c_str());
                    }
                    else {
                        if (ko_fluxes[i] <0)
                            pl_fillcolorname(colors[1].c_str());
                        else
                            pl_fillcolorname(colors[3].c_str());
                    }
                }
                
                if ((flux_vals[i] !=0) &&(ko_fluxes[i] !=0)){
                    
                    val =((log(fabs(flux_vals[i]))/log(10))-intlow)/(intlog-intlow);
                    val2=((log(fabs(ko_fluxes[i]))/log(10))-intlow)/(intlog-intlow);
                    
                   // cout<<"Log of flux for "<<rxn_names[i]<<" is "<<log(fabs(flux_vals[i]))/log(10)<<" and "<<(log(fabs(ko_fluxes[i]))/log(10))<<" V1 "<<val<<" V@: "<<val2<<" Flux "<<flux_vals[i]<<" KF:" <<ko_fluxes[i]<<endl;
                    
                    pl_fcircle (val*xrange+border+line_off, val2*yrange+border+line_off, dot_size);
                }
                else {
                    if (flux_vals[i] ==0) {
                        val2=((log(fabs(ko_fluxes[i]))/log(10))-intlow)/(intlog-intlow);
                        pl_fcircle (border, val2*yrange+border+line_off, dot_size);
                    }
                    if (ko_fluxes[i] ==0) {
                        val =((log(fabs(flux_vals[i]))/log(10))-intlow)/(intlog-intlow);
                        
                        pl_fcircle (val*xrange+border+line_off, border, dot_size);
                    }
                }
            }
            
            for(i=0; i<num_fluxes; i++) {
                //pl_pencolorname ("white");
                if (rxn_names[i]==biomass_name) {
                    pl_fillcolorname(colors[4].c_str());
                    
                    val =((log(fabs(flux_vals[i]))/log(10))-intlow)/(intlog-intlow);
                    val2=((log(fabs(ko_fluxes[i]))/log(10))-intlow)/(intlog-intlow);
                    pl_fcircle (val*xrange+border+line_off, val2*yrange+border+line_off, dot_size);
                }
            }
        }
        pl_pencolorname ("black");
        
        
    }
    else {
        pl_flinewidth (0.5);
        pl_pencolorname ("black");
        pl_fline(border+line_off, border+line_off, spacex-spacer, border+line_off);
        pl_fline(border+line_off, spacey-spacer, spacex-spacer, spacey-spacer);
        pl_fline(spacex-spacer, border+line_off, spacex-spacer, spacey-spacer);
        pl_fline(border+line_off, border+line_off, border+line_off, spacey-spacer);
        
        maxcnt=flux_cnts[0];
        for(i=1; i<11; i++) {
            if (flux_cnts[i]>maxcnt) maxcnt=flux_cnts[i];
        }
        //cout<<"Max # rxns: "<<maxcnt<<endl;
        newval=border+line_off+ (0.5/11.0)*xrange;
        
        pl_fmove(newval, 0.88*border);
        mark=(double)i;
        
        pl_alabel('c', 'c', "Lethal");
        pl_fline(newval, border+line_off, newval, 1.05*(border+line_off));
        
        
        for(i=1; i<=10; i++) {
            newval=((((double)i+0.5)/11.0)*xrange)+(border+line_off);
            
            pl_fmove(newval, 0.88*border);
            mark=(double)i;
            std::stringstream label_ss2;
            
            label_ss2<<((double)i-0.5)/10.0;
            pl_alabel('c', 'c', label_ss2.str().c_str());
            pl_fline(newval, border+line_off, newval, 1.05*(border+line_off));
        }
        
        logval=log(maxcnt)/log(10.0);
        intlog=(int)logval;
        
        y_topl=pow(10, intlog);
        double myinc;
        if (maxcnt > (0.5*pow(10, intlog+1))) myinc=0.5*pow(10, intlog);
        else myinc=0.25*pow(10, intlog);
        while(y_topl<maxcnt)
            y_topl+=myinc;
        
        if (y_topl == (double)maxcnt) y_topl+=myinc;
      
        ymark_l=y_topl/(double)ymarks;
        
        //cout<<"Top val: "<<y_topl<<" int log "<<intlog<<" # y: "<<ymarks<<" L: "<<ymark_l<<endl;
    
        
        for(i=0; i<=ymarks; i++) {
            // cout<<"Printing "<<label<<" at "<<newval<<" i "<<i<<" m "<<ymark_l<<" mark v "<<mark<<endl;
            newval=((double)i/ymarks)*(yrange)+(border+line_off);
            
            
            mark=(double)(i/ymarks);
            std::stringstream label_ss2;
            
            label_ss2<<ymark_l*i;
            
            pl_fmove(0.88*border-0.35*(pl_flabelwidth (label_ss2.str().c_str())), newval);
            pl_alabel('c', 'c', label_ss2.str().c_str());
            pl_fline(border+line_off, newval, 1.05*(border+line_off), newval);
        }
        pl_filltype(1);
        pl_fillcolorname(colors[1].c_str());
        
        for(i=0; i<11; i++) {
            //cout<<"Box "<<i<<"= "<<flux_cnts[i]<<endl;
            pl_fbox (border+(((double)i+0.05)/11.0)*xrange, border, border+((double)(i+0.95)/11.0)*xrange, border+(((double)flux_cnts[i])/y_topl)*yrange);
           
        }
        
    }

    if (pl_closepl () < 0)      /* close Plotter */ {
        fprintf (stderr, "Couldn't close Plotter\n");
        return 1;
      }

      pl_selectpl (0);            /* select default Plotter */
      if (pl_deletepl (thandle) < 0) /* delete Plotter we used */ {
          fprintf (stderr, "Couldn't delete Plotter\n");
          return 1;
      }
      
      if (streamplot==true) {
          //std::cout<<"Extracting from buffer\n";
          fflush(outfile);
          
          //printf ("At close buf = `%s', size = %zu\n", buf_pt, size);
          mypos=(size_t)0;
          
         
          while(mypos < size) {
              plot_ss->put((char)buf_pt[mypos]);
              //cout<<"P: "<<mypos<<" = "<<(int)buf_pt[mypos]<<endl;
              ++mypos;
          }
          
      }
 
    fclose(outfile);
    delete[] my_names;
    return 0;
}

