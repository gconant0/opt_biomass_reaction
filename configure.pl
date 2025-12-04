#!/usr/bin/perl

use strict;

my($i, $plot_path, $parallel, %paths, $plotlib_installed, $version, $cc, $CC, $make_daemon, $lapack_installed, $lapack_version, $blas_version, $plot_version, $cgi_installed, $boost_installed, @out, $myloc, $j);

$lapack_installed=0;
$boost_installed=0;
$plot_path ="";

$plotlib_installed=0;
$plot_version="";

if ($ARGV[0] =~ /plot/) {
    
    if((glob("/usr/lib/libplot*")) ) {
        $plotlib_installed=1;

        @out=`ls /usr/lib/libplot*`;
        for($i=0; $i<@out; $i++) {
            if ($out[$i] =~ /(libplot\.so\.\d)/) {
                $plot_version=$1;
            }
        }
    }
    
    if((glob("/usr/local/lib/libplot*")) ) {
        $plotlib_installed=1;

        @out=`ls /usr/local/lib/libplot*`;
        for($i=0; $i<@out; $i++) {
            if ($out[$i] =~ /(libplot\.so\.\d)/) {
                $plot_version=$1;
            }
        }
    }

    if((glob("/lib/libplot*")) ) {
        $plotlib_installed=1;

        @out=`ls /lib/libplot*`;
        for($i=0; $i<@out; $i++) {
            if ($out[$i] =~ /(libplot\.so\.\d)/) {
                $plot_version=$1;
            }
        }
    }
    
    if((glob("/usr/lib/x86_64-linux-gnu/libplot*")) ) {
        $plotlib_installed=1;

        @out=`ls /usr/lib/x86_64-linux-gnu/libplot*`;
        for($i=0; $i<@out; $i++) {
            if ($out[$i] =~ /(libplot\.so\.\d)/) {
                $plot_version=$1;
            }
        }
    }

    
}


if (-e "src/Makefile") {`rm src/Makefile`}

$lapack_version="";
$blas_version="";

if((glob("/usr/lib/liblapack*")) && (glob("/usr/lib/libblas*"))) {
    $lapack_installed=1;
    @out=`ls /usr/lib/liblapack*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(liblapack\.so\.\d)/) {
            $lapack_version=$1;
        }
    }
    
    @out=`ls /usr/lib/libblas*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(libblas\.so\.\d)/) {
            $blas_version=$1;
        }
    }
}

if((glob("/usr/local/lib/liblapack*")) && (glob("/usr/local/lib/libblas*"))) {
    $lapack_installed=1;
    
    @out=`ls /usr/local/lib/liblapack*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(liblapack\.so\.\d)/) {
            $lapack_version=$1;
        }
    }
    
    @out=`ls /usr/local/lib/libblas*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(libblas\.so\.\d)/) {
            $blas_version=$1;
        }
    }
    
}

if((glob("/lib/liblapack*")) && (glob("/lib/libblas*"))) {
    $lapack_installed=1;
    
    @out=`ls /lib/liblapack*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(liblapack\.so\.\d)/) {
            $lapack_version=$1;
        }
    }
    
    @out=`ls /lib/libblas*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(libblas\.so\.\d)/) {
            $blas_version=$1;
        }
    }
}

if((glob("/usr/lib/x86_64-linux-gnu/liblapack*")) && (glob("/usr/lib/x86_64-linux-gnu/libblas*"))) {
    $lapack_installed=1;
    
    @out=`ls /usr/lib/x86_64-linux-gnu/liblapack*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(liblapack\.so\.\d)/) {
            $lapack_version=$1;
        }
    }
    
    @out=`ls /usr/lib/x86_64-linux-gnu/libblas*`;
    for($i=0; $i<@out; $i++) {
        if ($out[$i] =~ /(libblas\.so\.\d)/) {
            $blas_version=$1;
        }
    }
}

$CC=`which g++`;
$cc=`which gcc`;

$CC =~ s/\n//;

if ($CC eq "") {
    $CC=`which c++`;
    if ($CC eq "") {
        $CC=`which xlC`;
        if ($CC eq "") {
            $CC=`which Clang`;
                if($CC eq "") {
                    $CC = `which icpc`;
                
                    if ($CC eq "") {
                        print "ERROR: No c++ compiler found\n";
                        exit();
                    }
                }
        }
	}
}

print "Found $CC for c++\n";

$cc=`which gcc`;

$cc =~ s/\n//;

if ($cc eq "") {
    $cc=`which cc`;
    if ($cc eq "") {
        $cc=`which xlc`;
        if ($cc eq "") {
                $cc=`which clang`;
                        if($cc eq "") {
                                $cc = `which icc`;

                                if ($cc eq "") {
                                        print "ERROR: No c compiler found\n";
                                        exit();
                                }
                        }
        }

    }
}

print "Found $cc for c\n";

if($lapack_installed == 0) {
	print "ERROR: Please install LAPACK first\n";
	exit();
}

open(WRITEMAKE, ">src/Makefile") or die;

print WRITEMAKE "#Makefile for opt_biomass_rxn\n";
print WRITEMAKE "#G. Conant, 11/29/25\n\n";
print WRITEMAKE "cc = $cc\n";
print WRITEMAKE "CC = $CC\n";

print WRITEMAKE "O = o\nSRC_DIR = ./\n";

print WRITEMAKE "CFLAGS = -DGCC_COMPILE \n";
print WRITEMAKE "INCLUDE = -I.\n";
print WRITEMAKE "OPTIM_SPEED = -O3\nOPTIM_SIZE = -O1\nMATH_LIB = -lm\n";
if ($lapack_version eq "") {
    print WRITEMAKE "LAPACK_LIB = -llapack\nBLAS_LIB = -lblas\n";
}
else {
    print WRITEMAKE "LAPACK_LIB = -l:", $lapack_version, "\nBLAS_LIB = -l:", $blas_version, "\n";
}

if ($plotlib_installed ==1) {
    if ($plot_version eq "") {
        print WRITEMAKE "PLOT_LIB = -l:", $plot_version, "\n";
    }
    else {
        print WRITEMAKE "PLOT_LIB = -lplot\n";
    }
}

print WRITEMAKE "OPTIONS = \$(CFLAGS) \$(INCLUDE)\n";
print WRITEMAKE "GLPK_LIB = -lglpk\n";

print WRITEMAKE "all:  opt_biomass_rxn\n";

print WRITEMAKE "OPT_BIOMASS_RXN_OBJS = opt_biomass_rxn.\$(O) stoich_mat.\$(O)  lin_program.\$(O)  gen_dna_funcs.\$(O) ";
if ($plotlib_installed ==1) {
    print WRITEMAKE " plot_fba.\$(O) ";
}
print WRITEMAKE "\n";


print WRITEMAKE "opt_biomass_rxn: \$(OPT_BIOMASS_RXN_OBJS)\n";
print WRITEMAKE "\t\$(CC) \$(LINUX_BUILD)  \$(LIBRARY_DIR)  -o ../opt_biomass_rxn \$(OPTIONS) \$(OPT_BIOMASS_RXN_OBJS)  \$(MATH_LIB)   \$(GLPK_LIB)  \$(LAPACK_LIB) \$(BLAS_LIB)";
if ($plotlib_installed ==1) {
    print WRITEMAKE " \$(PLOT_LIB) ";
}
print WRITEMAKE "\n";


print WRITEMAKE "%.o: %.cpp\n";
print WRITEMAKE "\t\$(CC)  \$(OPTIONS) \$(OPTIM_SPEED) -c \$<\n";


print WRITEMAKE "%.o: %.c\n";
print WRITEMAKE "\t\$(CC)  \$(OPTIONS) \$(OPTIM_SPEED) -c \$<\n";
close(WRITEMAKE);

open(WRITEMAKE, ">Makefile") or die;
print WRITEMAKE "#Master Makefile for opt_biomass_rx";


print WRITEMAKE "default: all\n";


    print WRITEMAKE "all:  \n";
    print WRITEMAKE "\tcd src;  make\n";











close(WRITEMAKE);
