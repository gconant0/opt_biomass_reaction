# opt_biomass_reaction
This package builds a tool for running flux-balance analysis, given a set of input reactions, a boundary file for nutrients and the name of the biomass reaction to maximize flux through. You can also specific reactions to remove prior to running the FBA analysis, or conduct a sequential removal of all reactions that show flux in the wild-type model.

Usage: opt_biomass_rxn <matrix file> <boundary file> Biomass_Rxn_Name ((constrain_Rxn) (constraint val)) (-k:RxnName)   (-sko) <filename>

If -sko is specified, all single knockouts are made and if "filename" is given, the results of this analysis are stored to that file.

The code depends on the GNU Linear programming toolkit and LAPACK: https://www.gnu.org/software/glpk/.
