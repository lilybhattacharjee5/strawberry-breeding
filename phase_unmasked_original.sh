 beagle_folder_root="/Users/lilybhattacharjee/Documents/beagle";
 inputs_root="shuffle_nans_inputs"; 
 outputs_root="shuffle_nans_outputs"; 
 original_folder="all_masked";
 for run in {1..1} ; do 
    d=(${beagle_folder_root}/${inputs_root}/${original_folder}/${run}_inputs/) ; 
    gunzip -r ${d} ;
    pattern="${d}/genotype.vcf" ; 
    files=( $pattern ) ; 
    tr -d " " < ${files[0]} > ${d}/cleaned_genotypes.vcf ; 
    echo ${d}/cleaned_genotypes.vcf; 
    mkdir ${beagle_folder_root}/${outputs_root}/${original_folder}/${run}_outputs/; 
    java -Xmx1g -jar beagle.28Jun21.220.jar gt=$d/cleaned_genotypes.vcf out=${beagle_folder_root}/${outputs_root}/${original_folder}/${run}_outputs/phased_genotypes seed=0 ;
    gunzip -r ${beagle_folder_root}/${outputs_root}/${original_folder}/${run}_outputs/; 
done 
