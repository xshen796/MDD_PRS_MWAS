module load igmm/apps/R/3.6.1

rm XS.log
rm -r /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/*

# Prepare exposure data  ----------------------------------------------------------------------

Rscript PREP.exposure_multivar.R --SNPlist ls.snp.multivar_MR.txt \
--CpGlist ls.cpg.multivar_MR.txt \
--mqtl /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/mQTL/mQTL_forMR_GS/ \
--interexp /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/exposure_stats/



# Run MR   ------------------------------------------------------------------------------------
    targetfile="/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/summstats/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.forPRSice.gz"
    result_tag=`echo $targetfile | sed "s/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.forPRSice/MDD/"`
    result_tag=`echo $result_tag | sed "s/\\/exports\\/igmm\\/eddie\\/GenScotDepression\\/shen\\/bakup.dat\\/summstats\\///"`  
    result_tag=`echo $result_tag | sed "s/.gz//"`   
    
    rm -r /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/outcome_stats/

    # Prepare outcome data 
    Rscript PREP.outcome_multivar.R \
    --exposure /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/exposure_stats/ \
    --summout ${targetfile} \
    --outcomeName ${result_tag} \
    --interout /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/outcome_stats/

    # MR analysis 
    Rscript ANALY.multivar_MR.R \
    --MRexposure /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/exposure_stats/multivar.exposure_dat \
    --MRoutcome /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/outcome_stats/multivar.outcome_dat \
    --outcomeName ${result_tag} \
    --outputfig /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/GS_MR_DNAm_to_MDD/DNAm_to_${result_tag}_mvMR_figs/ \
    --outputtable /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/GS_MR_DNAm_to_MDD/Summary_DNAm_to_${result_tag}_mvMR.txt \
    --saveHarmonisedData T
