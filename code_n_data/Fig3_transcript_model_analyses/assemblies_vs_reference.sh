#!/bin/bash
#PBS -l select=1:ncpus=10:mem=20G 

# load variables from general configuration file
if [ -z "$path" ]
then
      CURR_DIR=$1
else
      CURR_DIR=$path
fi					                                                # obtain current script directory

SCRIPTDIR="$CURR_DIR"
CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)                                    # obtain configuration file directory
source $CONFIG/general/config.sh
source $CONFIG/general/images.sh "$CONFIG/general/"

# load local configuration file
source $CURR_DIR/config.sh 

# directories
WD="$BASEDIR/analysis/benchmarking/assemblies_vs_reference"
GTF_LRGASP="$WD/50_genes_reference/LRGASP_final.gtf"
mkdir -p $WD/plots/tx_tables/gencode
mkdir -p $WD/plots/gffcompare_tables/LRGASP
mkdir -p $WD/plots/gffcompare_tables/gencode


### correct strandedness in unstranded files
#awk 'BEGIN {FS=OFS="\t"} {if ($7 == ".") $7="+";} 1' $WD/isoquant_sensitive/Neuron_Series_THP1.transcript_models.gtf > $WD/isoquant_sensitive/Neuron_Series_THP1.transcript_models.stranded_corrected.gtf
#awk 'BEGIN {FS=OFS="\t"} {if ($7 == ".") $7="+";} 1' $WD/isoquant_standard/Neuron_Series_THP1.transcript_models.gtf > $WD/isoquant_standard/Neuron_Series_THP1.transcript_models.stranded_corrected.gtf

### arrays for comparison with GENCODE genes
declare -a arr=("new_assembler" "talon" "isoquant_sensitive" "isoquant_standard" "bambu_NDR0169" "bambu_NDR02" "bambu_NDR05")
declare -a assembly=($FILES/RIKEN_files_January_2024/table2.noIP.detected.alone.5read.gtf $FILES/RIKEN_files_January_2024/camilla.TALON.table2.noIP.5read.gtf $WD/isoquant_sensitive/Neuron_Series_THP1.transcript_models.stranded_corrected.gtf $WD/isoquant_standard/Neuron_Series_THP1.transcript_models.stranded_corrected.gtf $WD/bambu_NDR0169/full_bambu_stranded_corrected_separated.gtf $WD/bambu_NDR02/full_bambu_stranded_corrected_separated.gtf $WD/bambu_NDR05/full_bambu_stranded_corrected_separated.gtf)
declare -a out_name=("new_assembler_vs_gencode" "talon_vs_gencode" "isoquant_permissive_vs_gencode" "isoquant_standard_vs_gencode" "bambu_NDR0169_vs_gencode" "bambu_NDR02_vs_gencode" "bambu_NDR05_vs_gencode")

### arrays for comparison with LRGASP genes
declare -a out_name_LRGASP=("new_assembler_vs_LRGASP" "talon_vs_LRGASP" "isoquant_permissive_vs_LRGASP" "isoquant_standard_vs_LRGASP" "bambu_NDR0169_vs_LRGASP" "bambu_NDR02_vs_LRGASP" "bambu_NDR05_vs_LRGASP")

### arrays for comparison with internal priming-cleaned assemblies
declare -a arr_IP=("SALA" "talon" "isoquant_sensitive" "isoquant_standard" "bambu_NDR0169" "bambu_NDR02" "bambu_NDR05")
files_IP=$(find $WD/internal_primed -type f)
assembly_IP_2=($files_IP)
assembly_IP=("$FILES/RIKEN_files_January_2024/table2.noIP.detected.alone.5read.gtf" "$FILES/RIKEN_files_January_2024/camilla.TALON.table2.noIP.5read.gtf")
assembly_IP+=("${assembly_IP_2[@]}")
declare -a out_name_IP=("SALA_IP_vs_gencode" "talon_IP_vs_gencode" "isoquant_permissive_IP_vs_gencode" "isoquant_standard_IP_vs_gencode" "bambu_NDR0169_IP_vs_gencode" "bambu_NDR02_IP_vs_gencode" "bambu_NDR05_IP_vs_gencode")

### arrays for comparison with internal priming-cleaned assemblies and SALA table 5
declare -a arr_IP=("SALA" "SALA_table5" "talon" "isoquant_sensitive" "isoquant_standard" "bambu_NDR0169" "bambu_NDR02" "bambu_NDR05")
files_IP=$(find $WD/internal_primed -type f)
assembly_IP_2=($files_IP)
assembly_IP=("$FILES/RIKEN_files_January_2024/table2.noIP.detected.alone.5read.gtf" "$SALA_TABLE5_ONLY_DETECTED" "$FILES/RIKEN_files_January_2024/camilla.TALON.table2.noIP.5read.gtf")
assembly_IP+=("${assembly_IP_2[@]}")
declare -a out_name_IP=("SALA_IP_vs_gencode" "SALA_table5_IP_vs_gencode" "talon_IP_vs_gencode" "isoquant_permissive_IP_vs_gencode" "isoquant_standard_IP_vs_gencode" "bambu_NDR0169_IP_vs_gencode" "bambu_NDR02_IP_vs_gencode" "bambu_NDR05_IP_vs_gencode")


### run comparison between GENCODE transcriptome and produced assemblies 
## now loop through the above array
#for index in "${!arr[@]}";do
	#echo "$index -> ${arr[$index]}"
	#mkdir -p $WD/"${arr[$index]}"/gffcompare/"${out_name[$index]}"
	#mkdir -p $WD/"${arr[$index]}"/SQANTI/"${out_name[$index]}"

	### Run gffcompare
	#$SINGC /gffcompare/gffcompare -o $WD/"${arr[$index]}"/gffcompare/"${out_name[$index]}"/"${out_name[$index]}" -r $GTF_ANNOT -R ${assembly[$index]}

	### Build tables for plots from gffcompare output
        #tail -n +11 $WD/"${arr[$index]}"/gffcompare/"${out_name[$index]}"/"${out_name[$index]}".stats| head -n 6 | awk -F'[:|]' '{gsub(/[ \t]+/, "", $1); gsub(/[ \t]+/, "", $2); gsub(/[ \t]+/, "", $3); print $1 "\t" $2 "\t" $3}' > $WD/plots/gffcompare_tables/gencode/"${arr[$index]}".txt

	### Run SQANTI3
	#source activate ${ENVS}/SQANTI3.env            # SQANTI3 was installed from the Conesa Lab github page
	#cd $WD/"${arr[$index]}"/SQANTI/"${out_name[$index]}"
	#python3 $SQANTI_DIR/sqanti3_qc.py -t $THREADS --report both --force_id_ignore --aligner_choice=minimap2 "${assembly[$index]}" $GTF_ANNOT $GENOME_FA -o "${out_name[$index]}"
	#conda deactivate

	### Build tables for plots from SQANTI output
	#tail -n +2 $WD/"${arr[$index]}"/SQANTI/"${out_name[$index]}"/"${out_name[$index]}"_classification.txt | cut -f6 | sort | uniq -c > $WD/plots/tx_tables/gencode/"${arr[$index]}"_classes_stats.txt
	#awk '$3=="transcript"' "${assembly[$index]}" | sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr[$index]}"_tx_number.txt
	#awk '$3=="gene"' "${assembly[$index]}" | sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr[$index]}"_gene_number.txt
	#awk '$3=="gene"' "${assembly[$index]}" | awk -F'gene_id "' '{print $2}' | awk -F'"' '{print $1}' | grep -v '^ENSG'| sort | uniq | wc -l  > $WD/plots/tx_tables/gencode/"${arr[$index]}"_gene_number_novel.txt
	#awk '$3=="gene"' "${assembly[$index]}" | awk -F'gene_id "' '{print $2}' | awk -F'"' '{print $1}' | grep '^ENSG'| sort | uniq | wc -l  > $WD/plots/tx_tables/gencode/"${arr[$index]}"_gene_number_annotated.txt
	#awk '$3=="transcript"' "${assembly[$index]}" | awk -F'transcript_id "' '{print $2}' | awk -F'"' '{print $1}' |  grep -v '^ENST'| sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr[$index]}"_tx_number_novel.txt
	#awk '$3=="transcript"' "${assembly[$index]}" | awk -F'transcript_id "' '{print $2}' | awk -F'"' '{print $1}' |  grep '^ENST'| sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr[$index]}"_tx_number_annotated.txt
#done


### run comparison between LRGASP consortium selected genes and produced assemblies
#for index in "${!arr[@]}";do
        #echo "$index -> ${arr[$index]}"
        #mkdir -p $WD/"${arr[$index]}"/gffcompare/"${out_name_LRGASP[$index]}"
        #mkdir -p $WD/"${arr[$index]}"/SQANTI/"${out_name_LRGASP[$index]}"

        ### Run gffcompare
        #$SINGC /gffcompare/gffcompare -o $WD/"${arr[$index]}"/gffcompare/"${out_name_LRGASP[$index]}"/"${out_name_LRGASP[$index]}" -r $GTF_LRGASP -R ${assembly[$index]}
	
	### Build tables for plots from gffcompare output
        #tail -n +11 $WD/"${arr[$index]}"/gffcompare/"${out_name_LRGASP[$index]}"/"${out_name_LRGASP[$index]}".stats| head -n 6 | awk -F'[:|]' '{gsub(/[ \t]+/, "", $1); gsub(/[ \t]+/, "", $2); gsub(/[ \t]+/, "", $3); print $1 "\t" $2 "\t" $3}' > $WD/plots/gffcompare_tables/LRGASP/"${arr[$index]}".txt

	### Run SQANTI3
        #source activate ${ENVS}/SQANTI3.env            # SQANTI3 was installed from the Conesa Lab github page
        #cd $WD/"${arr[$index]}"/SQANTI/"${out_name_LRGASP[$index]}"
        #python3 $SQANTI_DIR/sqanti3_qc.py -t $THREADS --report both --force_id_ignore --aligner_choice=minimap2 "${assembly[$index]}" $GTF_LRGASP $GENOME_FA -o "${out_name_LRGASP[$index]}"
        #conda deactivate

#done


### run comparison between GENCODE transcriptome and internal priming-cleaned assemblies
#for index in "${!arr_IP[@]}";do
#        echo "$index -> ${arr_IP[$index]}"
#	 echo "${assembly_IP[$index]}"
#        mkdir -p $WD/"${arr_IP[$index]}"/gffcompare/"${out_name_IP[$index]}"
#        mkdir -p $WD/"${arr_IP[$index]}"/SQANTI/"${out_name_IP[$index]}"

        ### Run gffcompare
        #$SINGC /gffcompare/gffcompare -o $WD/"${arr_IP[$index]}"/gffcompare/"${out_name_IP[$index]}"/"${out_name_IP[$index]}" -r $GTF_ANNOT -R ${assembly_IP[$index]}

        ### Build tables for plots from gffcompare output
        #tail -n +11 $WD/"${arr_IP[$index]}"/gffcompare/"${out_name_IP[$index]}"/"${out_name_IP[$index]}".stats| head -n 6 | awk -F'[:|]' '{gsub(/[ \t]+/, "", $1); gsub(/[ \t]+/, "", $2); gsub(/[ \t]+/, "", $3); print $1 "\t" $2 "\t" $3}' > $WD/plots/gffcompare_tables/gencode/"${arr_IP[$index]}".txt

        ### Run SQANTI3
        #source activate ${ENVS}/SQANTI3.env            # SQANTI3 was installed from the Conesa Lab github page
        #cd $WD/"${arr_IP[$index]}"/SQANTI/"${out_name_IP[$index]}"
        #python3 $SQANTI_DIR/sqanti3_qc.py -t $THREADS --report both --force_id_ignore --aligner_choice=minimap2 "${assembly_IP[$index]}" $GTF_ANNOT $GENOME_FA -o "${out_name_IP[$index]}"
        #conda deactivate

        ### Build tables for plots from SQANTI output
        #tail -n +2 $WD/"${arr_IP[$index]}"/SQANTI/"${out_name_IP[$index]}"/"${out_name_IP[$index]}"_classification.txt | cut -f6 | sort | uniq -c > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_classes_stats.txt
        #awk '$3=="transcript"' "${assembly_IP[$index]}" | sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_tx_number.txt
        #awk '$3=="gene"' "${assembly_IP[$index]}" | sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_gene_number.txt
        #awk '$3=="gene"' "${assembly_IP[$index]}" | awk -F'gene_id "' '{print $2}' | awk -F'"' '{print $1}' | grep -v '^ENSG'| sort | uniq | wc -l  > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_gene_number_novel.txt
        #awk '$3=="gene"' "${assembly_IP[$index]}" | awk -F'gene_id "' '{print $2}' | awk -F'"' '{print $1}' | grep '^ENSG'| sort | uniq | wc -l  > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_gene_number_annotated.txt
        #awk '$3=="transcript"' "${assembly_IP[$index]}" | awk -F'transcript_id "' '{print $2}' | awk -F'"' '{print $1}' |  grep -v '^ENST'| sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_tx_number_novel.txt
        #awk '$3=="transcript"' "${assembly_IP[$index]}" | awk -F'transcript_id "' '{print $2}' | awk -F'"' '{print $1}' |  grep '^ENST'| sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_tx_number_annotated.txt

#done



### run comparison between GENCODE transcriptome and internal priming-cleaned assemblies and SALA table5
for index in "${!arr_IP[@]}";do
        echo "$index -> ${arr_IP[$index]}"
        echo "${assembly_IP[$index]}"
        mkdir -p $WD/"${arr_IP[$index]}"/gffcompare/"${out_name_IP[$index]}"
        mkdir -p $WD/"${arr_IP[$index]}"/SQANTI/"${out_name_IP[$index]}"

        ### Run gffcompare
        #$SINGC /gffcompare/gffcompare -o $WD/"${arr_IP[$index]}"/gffcompare/"${out_name_IP[$index]}"/"${out_name_IP[$index]}" -r $GTF_ANNOT -R ${assembly_IP[$index]}

        ### Build tables for plots from gffcompare output
        #tail -n +11 $WD/"${arr_IP[$index]}"/gffcompare/"${out_name_IP[$index]}"/"${out_name_IP[$index]}".stats| head -n 6 | awk -F'[:|]' '{gsub(/[ \t]+/, "", $1); gsub(/[ \t]+/, "", $2); gsub(/[ \t]+/, "", $3); print $1 "\t" $2 "\t" $3}' > $WD/plots/gffcompare_tables/gencode/"${arr_IP[$index]}".txt

        ### Run SQANTI3
        source activate ${ENVS}/SQANTI3.env            # SQANTI3 was installed from the Conesa Lab github page
        cd $WD/"${arr_IP[$index]}"/SQANTI/"${out_name_IP[$index]}"
        python3 $SQANTI_DIR/sqanti3_qc.py -t $THREADS --report both --force_id_ignore --aligner_choice=minimap2 "${assembly_IP[$index]}" $GTF_ANNOT $GENOME_FA -o "${out_name_IP[$index]}"
        conda deactivate

        ### Build tables for plots from SQANTI output
        tail -n +2 $WD/"${arr_IP[$index]}"/SQANTI/"${out_name_IP[$index]}"/"${out_name_IP[$index]}"_classification.txt | cut -f6 | sort | uniq -c > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_classes_stats.txt
        awk '$3=="transcript"' "${assembly_IP[$index]}" | sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_tx_number.txt
        awk '$3=="gene"' "${assembly_IP[$index]}" | sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_gene_number.txt
        awk '$3=="gene"' "${assembly_IP[$index]}" | awk -F'gene_id "' '{print $2}' | awk -F'"' '{print $1}' | grep -v '^ENSG'| sort | uniq | wc -l  > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_gene_number_novel.txt
        awk '$3=="gene"' "${assembly_IP[$index]}" | awk -F'gene_id "' '{print $2}' | awk -F'"' '{print $1}' | grep '^ENSG'| sort | uniq | wc -l  > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_gene_number_annotated.txt
        awk '$3=="transcript"' "${assembly_IP[$index]}" | awk -F'transcript_id "' '{print $2}' | awk -F'"' '{print $1}' |  grep -v '^ENST'| sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_tx_number_novel.txt
        awk '$3=="transcript"' "${assembly_IP[$index]}" | awk -F'transcript_id "' '{print $2}' | awk -F'"' '{print $1}' |  grep '^ENST'| sort | uniq | wc -l > $WD/plots/tx_tables/gencode/"${arr_IP[$index]}"_tx_number_annotated.txt

done
