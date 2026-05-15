# Vaccinium ALMT comparative-genomics

Data, analysis scripts, and intermediate outputs supporting the manuscript:

**Compact retention and regulatory-region positive selection of ALMT genes in acidophilic Vaccinium**
Bin Li, Xianyang Zhao, Ruiyi Fan (2026), Plant Cell Reports (in submission)

## Contents

- `scripts/` - Figure-generation and analysis scripts
- `data/alignments/` - MAFFT protein alignment + trimAl-trimmed; codon alignments for PAML
- `data/trees/` - IQ-TREE consensus tree, ML tree, iTOL annotation files
- `data/orthogroups/` - OrthoFinder Orthogroups.txt + verified ALMT IDs
- `data/kaks/` - Pairwise Ka/Ks tables (per orthogroup + per species)
- `data/paml/` - PAML branch-site outputs for Subfamily 5 (OG0008930) and Subfamily 2 (OG0002701)
- `data/synteny/` - MCScanX collinearity files + GFF inputs
- `data/structure/` - AlphaFold2-predicted Vduc1 PDB model

## Software versions

See manuscript Methods (Section 2) for detailed software versions and parameters.
Briefly: MAFFT v7.526, trimAl, IQ-TREE 2.3.6, OrthoFinder v3.1.4,
KaKs_Calculator v2.0, PAML 4.10, MCScanX, ColabFold/AlphaFold2.

## Citing

If you use any data or script from this repository, please cite the manuscript.

## Contact

Corresponding author: Prof. Dr. Ruiyi Fan, fanruiyi@outlook.com
