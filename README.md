# Biogeographic dating of the Tetragnatha radiations on Hawaii 
This work is yet unpublished. For more questions, please get in touch!

`prepareGeogRangeData.R` 	: Produces a Nexus file from biogeographic data
`partitionAlignment.py` 	: Produces separate locus alignments from a concatenated a;ignment
`nonclock.Rev`				: Performs a non-clock analysis for use as a starting tree
`biogeogDating.Rev` 		: Joint biogeographical and relaxed clock model for dating


## Constraints
Root age ~ U(22.8, 34.0)

## Notes:
- Landis' code produces 3 types of files
    - 1) parameter estimates for substitution model, relaxed clock and birth death models "run.log"
    - 2) stochastic character maps ("stoch.log")
    - 3) ancestral range maps ("states.log")

- Stochastic character mapping by rejection sampling (i.e, evolve states forward from root, only keep simulations with the same observed trait at tips; Nielsen 2002)


## Files from Landis
`make_anc_state.Rev`
- Imports tree trace, calculates the *maximum a posteriori* (MAP) tree, the consensus tree and the maximi, clade credibility tree.
- Imports state log using 


## To do
- Inter-island distances (what is the distance used by landis to continental areas?)
- We will probably have to do model testing
	- Distance-weighted vs. not
	- UCLN clock vs. uexp?


Composes of the following models

## Model 1 (M1):
* No biogeographic model + uniform root age prior of 0.1 - 15 Mya (split between viridis+versicolor and Hawaiians)

## Model 2 (M2)
** Biogeographic model + uniform root age prior of 0.1 - 15 Mya (split between viridis+versicolor and Hawaiians)

## Model 3 (M3)
** Clock model + no biogeography

## Model 4 (M4)
** Clock model + biogeographic model

## Loci used
Mitochondrial
- COI (2905-3322) coding
- CytB (2184-2539) coding
- 16S (2540-2904) coding
- 12S (1166-1524) coding

Nuclear
- ITS (1853-2183) non-coding
- 18S (1-269)
- SSU (796-1165)
2. Actin (270-495)
3. 28S (496-795)
6. H3 (1525-1852)