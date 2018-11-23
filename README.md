# Biogeographic dating of the Tetragnatha radiations on Hawaii 
This work is yet unpublished. For more questions, please get in touch!


`prepareGeogRangeData.R` 	: Produces a Nexus file from biogeographic data
`partitionAlignment.py` 	: Produces separate locus alignments from a concatenated a;ignment
`nonclock.Rev`				: Performs a non-clock analysis for use as a starting tree
`biogeogDating.Rev` 		: Joint biogeographical and relaxed clock model for dating


# Constraints
Root age ~ U(22.8, 34.0)

# Notes:
- Landis' code produces 3 types of files
    - 1) parameter estimates for substitution model, relaxed clock and birth death models "run.log"
    - 2) stochastic character maps ("stoch.log")
    - 3) ancestral range maps ("states.log")

- Stochastic character mapping by rejection sampling (i.e, evolve states forward from root, only keep simulations with the same observed trait at tips; Nielsen 2002)



# To do
- Inter-island distances (what is the distance used by landis to continental areas?)
- We will probably have to do model testing
	- Distance-weighted vs. not
	- UCLN clock vs. uexp?