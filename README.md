# snp_call_by_maping
maping reads to reference and call snp  
Note that snp position/site caused by snv site homologous recombination not removed.

## Important arguments  

* -ref_genome: the path of reference genome  
* -target_reads: the directory path which hold all reads file of all strains. the reads file belong to same strain should put into a sub-directory.  
* -snp_coverage_cutoff: cutoff value to filter snp position which contain strain whose base was lowly covered.   
* -snp_qual_cutoff: cutoff value to filter snp position which contain strain whose snp site have low QUAL value.  
* -snp_diversity_cutoff: cutoff value to filter low diversified snp position.

## Strategy
1. For reads files of each strain using bwa maping reads date to reference genome.  
2. Use samtools and bcftools to mpileup and calling snp of each strain.  
3. Union snp site(base is different to reference base) information of each strain and get all candidate snp position(position according to reference genome which have one or more varition site find in item strains). Note that every snp position contain many bases each base belong to a strain. And some base different to reference base.
4. Filter snp position using snp_coverage_cutoff, snp_qual_cutoff and snp_diversity_cutoff.  Filter method described as follow:
>* For each snp position, if there is one or more strain whose base's coverage is below snp_coverage_cutoff value this snp position will be removed.
>* For each snp position, if there is one or more strain whose base's (this is a snp site under context here) is less than snp_qual_cutoff value this snp position will be removed.
>* For each snp position, if diversity value is less than snp_diversity_cutoff, which will be removed. the snp diversity was defined as: 1-(most_present_base_amount/strain_amount).
5. After filter procedure all snp information will be printed with fasta format. which contain snp position, snp refrence base, base of each strain.

## Help information
`call_snp.py -h` will print help information.
