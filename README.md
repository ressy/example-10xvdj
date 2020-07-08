# 10x Example - VDJ Analysis

## Setup

    conda env update --file environment.yml
    conda activate example-10xvdj

(Requires [conda](https://www.anaconda.com/).)

Download Cell Ranger from
<https://support.10xgenomics.com/single-cell-vdj/software/downloads/latest?>
(can't give a direct link as it includes an access token in the URL) and put
the extracted directory in `$PATH`.

Then:

    snakemake -j

[10x.csv](data/10x.csv) links to the files listed here:
<https://support.10xgenomics.com/single-cell-vdj/datasets/3.1.0/vdj_v1_hs_pbmc3>

> PBMCs of a Healthy Donor (v1)
> Single Cell Immune Profiling Dataset by Cell Ranger 3.1.0

...as well as one GRCh38 reference.  Snakemake will download the ones it needs.

## Journey of one read

For example, `A00519:312:HCVFYDRXX:2:1151:32714:23547`.

### Raw Read Data

From the lane 2 files in `analysis/vdj_v1_hs_pbmc3_b_26x91_fastqs/`:

R1:

    @A00519:312:HCVFYDRXX:2:1151:32714:23547 1:N:0:GTCCGGTC
    AACCATGCAGTGACAGCACTATTCGC
    +
    ,:FF,FFF:FFF:F:FF,FF:,,FFF

R2:

    @A00519:312:HCVFYDRXX:2:1151:32714:23547 2:N:0:GTCCGGTC
    GCCCGAGTAGCAGGAGGAAGAGAAGCTGCGCGGGGGCTTCCATGGTTCCGTCTGGGTCCTAACTGAGCAGTTCCTCCCCAGCGAAGAAAGC
    +
    FFF,F:F:,FFFFFFFFF,FFFF,FF:FFFF,FFFFF,,FFF:FFF:FFFFF:FFF:FFFF:FFFFFFFF,FFFFFFFF:,,:,:FF:FFF

I1:

    @A00519:312:HCVFYDRXX:2:1151:32714:23547 1:N:0:GTCCGGTC
    GTCCGGTC
    +
    FFFFFF,F

All together:

    Cell barcode from first part of R1: AACCATGCAGTGACAG
    UMI from second part of R1: CACTATTCGC
    Sequence from reverse-complement of R2: GCTTTCTTCGCTGGGGAGGAACTGCTCAGTTAGGACCCAGACGGAACCATGGAAGCCCCCGCGCAGCTTCTCTTCCTCCTGCTACTCGGGC
    Index barcode from I1: GTCCGGTC

The UMI links this read to one original cDNA.  The cell barcode links this cDNA
to one cell.  I think the index barcode is used for Illumina sample
demultiplexing but in the example files this is already done.

### Aligned to assembled contig for this cell

Now from the cellranger output in `analysis/vdj_v1_hs_pbmc3_b_26x91/outs/`:

The sequence shows up in `all_contig.bam` aligned to the assembled contig
`AACCATGCAGTGACAG-1_contig_1`, along with the other reads under this contig.
Manually re-creating with the above:

    >AACCATGCAGTGACAG-1_contig_1
    ----------CTGGGGAGGAACTGCTCAGTTAGGACCCAGACGGAACCATGGAAGCCCCAGCGCAGCTTCTCTTCCTCCTGCTACTCTGGCTCCCAGATAC...
    >A00519:312:HCVFYDRXX:2:1151:32714:23547_amplicon
    GCTTTCTTCGCTGGGGAGGAACTGCTCAGTTAGGACCCAGACGGAACCATGGAAGCCCCCGCGCAGCTTCTCTTCCTCCTGCTACTCGGGC-------------

Or see [igv.png](igv.png) with this read shown on the upper-left and in the
extra window.

There are 302 reads in the BAM assigned to this contig for this cell barcode,
under two different UMIs (so, two physically distinct but equivalent cDNA
molecules to start with, leading to a few hundred sequenced amplicons).
These are visible in the optional fields in the BAM (`UB:Z:CACTATTCGC UR:Z:CACTATTCGC`
for this one and the other amplicons starting from that cDNA, and
`UB:Z:GCCAGGGTTA UR:Z:GCCAGGGTTA` for the second; I don't know why there
are two fields for it).

### Aggregated with other cells under one clonotype

`AACCATGCAGTGACAG-1_contig_1` is the only contig identified for this cell
(there could be more than one) and is one of seven cells associated with
`clonotype1`.

`all_contig_annotations.csv` shows the seven contigs under `clonotype1` for the
seven different cell barcodes, all grouped under a single consensus for this
clonotype.  They are all identified as chain IGK.

    barcode            is_cell contig_id                   high_confidence length chain v_gene   d_gene j_gene c_gene full_length productive cdr3         cdr3_nt                              reads umis raw_clonotype_id raw_consensus_id       
    AACCATGCAGTGACAG-1 True    AACCATGCAGTGACAG-1_contig_1 True            559    IGK   IGKV3-15 None   IGKJ2  IGKC   True        True       CQQYDNWPPYTF TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT 302   2    clonotype1       clonotype1_consensus_1 
    CAGCTAACACCGATAT-1 True    CAGCTAACACCGATAT-1_contig_1 True            549    IGK   IGKV3-15 None   IGKJ2  IGKC   True        True       CQQYDNWPPYTF TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT 274   2    clonotype1       clonotype1_consensus_1 
    CCATGTCAGACCGGAT-1 True    CCATGTCAGACCGGAT-1_contig_1 True            556    IGK   IGKV3-15 None   IGKJ2  IGKC   True        True       CQQYDNWPPYTF TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT 324   2    clonotype1       clonotype1_consensus_1 
    CGTTCTGAGCTCAACT-1 True    CGTTCTGAGCTCAACT-1_contig_1 True            565    IGK   IGKV3-15 None   IGKJ2  IGKC   True        True       CQQYDNWPPYTF TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT 187   3    clonotype1       clonotype1_consensus_1 
    CTGAAGTGTGAAGGCT-1 True    CTGAAGTGTGAAGGCT-1_contig_1 True            554    IGK   IGKV3-15 None   IGKJ2  IGKC   True        True       CQQYDNWPPYTF TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT 250   3    clonotype1       clonotype1_consensus_1 
    CTTTGCGCAAGGACTG-1 True    CTTTGCGCAAGGACTG-1_contig_1 True            479    IGK   IGKV3-15 None   IGKJ2  IGKC   True        True       CQQYDNWPPYTF TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT 197   2    clonotype1       clonotype1_consensus_1 
    TATTACCGTCGCTTCT-1 True    TATTACCGTCGCTTCT-1_contig_1 True            565    IGK   IGKV3-15 None   IGKJ2  IGKC   True        True       CQQYDNWPPYTF TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT 299   2    clonotype1       clonotype1_consensus_1

`consensus_annotations.csv` has a row for this one clonotype consensus, reporting the 1833 reads total coming from 16 UMIs.

    clonotype_id consensus_id           length chain v_gene   d_gene j_gene c_gene full_length productive cdr3         cdr3_nt                              reads umis 
    clonotype1   clonotype1_consensus_1 565    IGK   IGKV3-15 None   IGKJ2  IGKC   True        True       CQQYDNWPPYTF TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT 1833  16 

`consensus_annotations.json` gives more detail in JSON format including the seven cells belonging to this clonotype.

        {
            "aa_sequence": "MEAPAQLLFLLLLWLPDTTGEIVMTQSPATLSVSPGERATLSCRASQSVSSNLAWYQQKPGQAPRLLIYGTSTRATGIPARFSGSGSGTEFTLTISSLQSEDFAVYYCQQYDNWPPYTFGQGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDN", 
            "annotations": [
                {
                    "annotation_length": 97, 
                    "annotation_match_end": 97, 
                    "annotation_match_start": 53, 
                    "cigar": "44M521S", 
                    "contig_match_end": 44, 
                    "contig_match_start": 0, 
                    "feature": {
                        "chain": "IGK", 
                        "display_name": "IGKV3-15", 
                        "feature_id": 283, 
                        "gene_name": "IGKV3-15", 
                        "region_type": "5'UTR"
                    }, 
                    "mismatches": [], 
                    "score": 88
                }, 
                {
                    "annotation_length": 345, 
                    "annotation_match_end": 345, 
                    "annotation_match_start": 0, 
                    "cigar": "44S345M176S", 
                    "contig_match_end": 389, 
                    "contig_match_start": 44, 
                    "feature": {
                        "chain": "IGK", 
                        "display_name": "IGKV3-15", 
                        "feature_id": 284, 
                        "gene_name": "IGKV3-15", 
                        "region_type": "L-REGION+V-REGION"
                    }, 
                    "mismatches": [], 
                    "score": 675
                }, 
                {
                    "annotation_length": 39, 
                    "annotation_match_end": 39, 
                    "annotation_match_start": 0, 
                    "cigar": "390S39M136S", 
                    "contig_match_end": 429, 
                    "contig_match_start": 390, 
                    "feature": {
                        "chain": "IGK", 
                        "display_name": "IGKJ2", 
                        "feature_id": 215, 
                        "gene_name": "IGKJ2", 
                        "region_type": "J-REGION"
                    }, 
                    "mismatches": [], 
                    "score": 63
                }, 
                {
                    "annotation_length": 320, 
                    "annotation_match_end": 136, 
                    "annotation_match_start": 0, 
                    "cigar": "429S136M", 
                    "contig_match_end": 565, 
                    "contig_match_start": 429, 
                    "feature": {
                        "chain": "IGK", 
                        "display_name": "IGKC", 
                        "feature_id": 213, 
                        "gene_name": "IGKC", 
                        "region_type": "C-REGION"
                    }, 
                    "mismatches": [], 
                    "score": 272
                }
            ], 
            "barcode": null, 
            "cdr3": "CQQYDNWPPYTF", 
            "cdr3_seq": "TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT", 
            "cdr3_start": 365, 
            "cdr3_stop": 401, 
            "clonotype": "clonotype1", 
            "contig_name": "clonotype1_consensus_1", 
            "filtered": true, 
            "frame": null, 
            "high_confidence": true, 
            "info": {
                "cell_contigs": [
                    "CCATGTCAGACCGGAT-1_contig_1", 
                    "TATTACCGTCGCTTCT-1_contig_1", 
                    "AACCATGCAGTGACAG-1_contig_1", 
                    "CTTTGCGCAAGGACTG-1_contig_1", 
                    "CGTTCTGAGCTCAACT-1_contig_1", 
                    "CAGCTAACACCGATAT-1_contig_1", 
                    "CTGAAGTGTGAAGGCT-1_contig_1"
                ], 
                "cells": [
                    "AACCATGCAGTGACAG-1", 
                    "CAGCTAACACCGATAT-1", 
                    "CCATGTCAGACCGGAT-1", 
                    "CGTTCTGAGCTCAACT-1", 
                    "CTGAAGTGTGAAGGCT-1", 
                    "CTTTGCGCAAGGACTG-1", 
                    "TATTACCGTCGCTTCT-1"
                ], 
                "clonotype_freq": 7, 
                "clonotype_prop": 0.009067357512953367
            }

### Final aggregation across antibody loci

`clonotypes.csv` shows clonotype information aggregated across multiple loci
where possible, giving linked heavy and light chains.

clonotype1's contigs only cover IGK, so we don't know what IGH looks like for
cells in this lineage.

    clonotype_id frequency proportion       cdr3s_aa         cdr3s_nt                                 
    clonotype1   7         0.00906735751295 IGK:CQQYDNWPPYTF IGK:TGTCAGCAGTATGATAACTGGCCTCCGTACACTTTT

clonotype6, for another example, shows information for both IGH and IGK.

    clonotype_id frequency proportion       cdr3s_aa                            cdr3s_nt                                                                                
    clonotype6   2         0.00259067357513 IGH:CARPGTTGTTGLKNW;IGK:CQQYNNWPLTF IGH:TGTGCGAGACCCGGTACAACTGGAACGACGGGTTTAAAAAACTGG;IGK:TGTCAGCAGTATAATAACTGGCCTCTCACCTTC

There are 760 clonotype entries total, and summing the frequency column gives the 772
cells reported in `web_summary.html`.

See also:

 * <https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/vdj>
