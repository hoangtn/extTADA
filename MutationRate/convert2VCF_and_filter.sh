chr=$1

outDir=$(pwd)
#########CONVERT TO VCF FILES
fileName=Reduced.All_coordinate_chr${chr}.txt.annotated.snp_fw_seq3.extractColumns

cat ${outDir}/${fileName}|awk '{print $1,toupper($2)}'|sed 's/chr//g'|\
	awk '{print $1":"substr($2,2,1)":"substr($2,6,1), $2}'|sort > ${outDir}/${fileName}.temp1.txt

fileName2=All_coordinate_chr${chr}.txt

cat ${outDir}/${fileName2} |awk '$4!=$3 {print $2":"$3":"$4,$5}'|\
	sort > ${outDir}/${fileName}.temp2.txt

## Join to a vcf file

join -1 1 -2 1 ${outDir}/${fileName}.temp1.txt ${fileName}.temp2.txt |uniq |\
	sed 's/:/ /g'|awk '{print $1,$2, ".",$3,$4,$6,$6,$5}'|\
	tr " " "\t" |sort -k 2 -n|uniq > ${outDir}/${fileName}.vcf

rm ${outDir}/${fileName}.temp*

###########
###ANNOTATE VCF FILES
softwareDir=/hpc/users/nguyet26/InstallSoftware
dbDir=/hpc/users/nguyet26/psychen/resources/dbNSFP/MutationFileBasedPrediction

${softwareDir}/java/jre1.8.0_65/bin/java -jar \
${softwareDir}/snpSift/snpEff/SnpSift.jar dbnsfp -v -db ${dbDir}/db.shortened.dbnsfp31a.gz $outDir/${fileName}.vcf \
	> $outDir/${fileName}.vcf.annotated

vcfFile=${fileName}.vcf.annotated
out="Extract.damaging"

snpsift=$(echo /hpc/users/nguyet26/InstallSoftware/java/jre1.8.0_65/bin/java -jar /hpc/users/nguyet26/InstallSoftware/snpSift/snpEff/SnpSift.jar)

$snpsift filter "(dbNSFP_SIFT_pred =~ 'D') & (dbNSFP_Polyphen2_HDIV_pred =~ 'D') & (dbNSFP_Polyphen2_HVAR_pred =~ 'D') & (dbNSFP_LRT_pred =~ 'D') & (dbNSFP_PROVEAN_pred =~ 'D')" ${outDir}/${vcfFile} \
	> ${outDir}/${out}.5methods.${vcfFile}

pred=MutationTaster
$snpsift filter "(dbNSFP_${pred}_pred =~ '[AD]')" ${outDir}/${out}.5methods.${vcfFile}  \
	> ${outDir}/${out}.6methods.${vcfFile}


pred=MutationAssessor
$snpsift filter "(dbNSFP_${pred}_pred =~ '[HM]')" ${outDir}/${out}.6methods.${vcfFile}  \
	> ${outDir}/${out}.7methods.${vcfFile}
