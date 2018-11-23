# Modified susan's original biogeography file to remove molokai and add a continent column
# Used R script to convert into nexus file, and then opened the file and manually change type to "Standard" from "DNA"

# 



For now:

# Constraints
Root age ~ U(22.8, 34.0)

# Notes:
- Landis' code produces 3 types of files
    - 1) parameter estimates for substitution model, relaxed clock and birth death models "run.log"
    - 2) stochastic character maps ("stoch.log")
    - 3) ancestral range maps ("states.log")


- Stochastic character mapping by rejection sampling (i.e, evolve states forward from root, only keep simulations with the same observed trait at tips; Nielsen 2002)



To do:
- Inter-island distances (what is the distance used by landis to continental areas?)
- We will probably have to do model testing
	- Distance-weighted vs. not
	- UCLN clock vs. uexp?