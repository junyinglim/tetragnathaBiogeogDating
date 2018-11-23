# partitionAlignmented.py
# Created: 17th October 2018
# Script splits the matK-trnK intron multiple sequence alignment by Wanke et al (2007) into coding and non-coding regions for partitioned phylogenetic analysis

import os
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

seq_dir = "/Users/junyinglim/Dropbox/spiders/phylo/ultros"
output_dir = "/Users/junyinglim/Dropbox/spiders/phylo/tetragnathaBiogeogDating/data/alignments"

full_align = AlignIO.read(os.path.join(seq_dir, "ultros.fas"), "fasta")
    
# Locus positions
# 1. 18S (1-269) AGCGCACACAAGAAA
# 2. Actin (270-495) TGCTGCATCTTC
# 3. 28S (496-795) CTGTGGGATGAA
# 4. SSU (796-1165) ATGCATGTCTAAG
# 5. 12S (1166-1524) ACCTTA
# 6. H3 (1525-1852) TCGTAAGAGTA
# 7. ITS (1853-2183) AAGAACGCAGCC
# 8. CytB (2184-2539) TCTTTTATC
# 9. 16S (2540-2904) AAAATTTAAAGG
# 10. COI (2905-3322) TATAAATAATTTAAG

locusNames = ["18S", "Actin", "28S", "SSU", "12S", "H3", "ITS", "CytB", "16S", "COI"]
locusStart = [0,269,495,795,1165,1524,1852,2183,2539,2904]
locusEnd = [269,495,795,1165,1524,1852,2183,2539,2904,3322]


for i in range(0, len(locusNames)):
	locusAlign = [seq[locusStart[i]:locusEnd[i]] for seq in full_align]
	msa = MultipleSeqAlignment(locusAlign)
	AlignIO.write(msa, os.path.join(output_dir, "ultros_"+locusNames[i]+".fasta"), "fasta")

