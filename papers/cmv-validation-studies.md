# Published CMV NGS Validation Studies

## Mallory et al. 2023 (PMCID: PMC10729743)
- **Title**: Development and validation of a next-generation sequencing assay with open-access analysis software for detecting resistance-associated mutations in CMV
- **Journal**: J Clin Microbiol, 2023;61(12):e00829-23
- **PMID**: 38092673
- **Assay**: Amplicon-based Ion Torrent NGS
- **Target genes**: UL27, UL54, UL56, UL97
- **Coverage**: 83 overlapping amplicons, ~10 kb, 138 mutation sites
- **Validation**: 
  - Custom plasmids with up to 99 known resistance mutations (Merlin backbone, GenBank NC_006273.2)
  - 185 patient specimens reanalyzed by NGS vs Sanger
  - NGS detected same mutations as Sanger in 68/69 (98.6%) specimens
  - 18 resistance mutations detected by NGS that Sanger missed (14 low-frequency <20%, 6 changed drug interpretation)
- **Database**: ARUPCMVResistanceCaller SQLite database (250 curated mutations)

## Streck et al. 2023 (PMCID: PMC10595055)
- **Title**: Use of next-generation sequencing to detect mutations associated with antiviral drug resistance in cytomegalovirus
- **Journal**: J Clin Microbiol, 2023;61(10):e00429-23
- **PMID**: 37750719
- **Assay**: NGS method for UL54/UL97 resistance mutations
- **Samples**: 27 plasma samples vs Sanger
- **Results**: 
  - 100% overall agreement for resistance/susceptibility determination (27/27)
  - 88% agreement at the level of resistance-associated mutations (22/25)
  - LOD: 500 IU/mL
  - Variant detection threshold: 15%
- **Software**: DeepChek-CMV (ABL Diagnostics)

## ARUPCMVResistanceCaller Database
- **Source**: GitHub mmallory26/ARUPCMVResistanceCaller, v2.1.0.5 release
- **Database**: CMV_resistance_mutations_20220321.db (SQLite)
- **Content**: 250 curated CMV resistance mutations
- **Columns**: gene, variant, codon_start, codon_end, ref_aa, var_aa, cidofovir, foscarnet, ganciclovir, maribavir, letermovir
- **Genes**: UL27, UL51, UL54, UL56, UL89, UL97
- **Drug classes**: 5 (ganciclovir, foscarnet, maribavir, letermovir, cidofovir)

## Benchmark Results
- 127/250 ARUP mutations overlap with HerpesDRG BED positions
- 100% vcfanno annotation match on all 127 benchmark mutations
- All 6 CMV resistance genes covered
- All 5 drug classes covered