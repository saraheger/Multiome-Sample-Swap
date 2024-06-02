# Multiome Sample Swap
Identity check 10x Multiome samples

Full article on detecting sample swaps:
https://gatk.broadinstitute.org/hc/en-us/articles/360041696232-Detecting-sample-swaps-with-Picard-tools

STEP 1: `ExtractFingerprint` (gatk) to extract fingerprint from a cellranger BAM
This makes the STEP 2 faster because the fingerprints are VCFs.

STEP 2: `CrosscheckFingerprints` across multiome GEX (still need to add ATAC)