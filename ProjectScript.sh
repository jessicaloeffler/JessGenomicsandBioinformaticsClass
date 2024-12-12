threads=8

#List SRR numbers
srr_files=("SRR11805799" "SRR11805801" "SRR11805803" "SRR11805805" "SRR11805807" "SRR11805809")

# Loop over each SRR file and run fasterq-dump with the specified number of threads
for srr in "${srr_files[@]}"; do
 echo "Downloading $srr with $threads threads..."
  fasterq-dump -t $threads $srr 
	done
	
echo "ALL DONE"

#after this, i went in an manually changed the names of the SRR files so that I didnt have to use the numbers. I'm sure there was a way to do it 
#in here, but I'm new and not that savvy. 
#SRR11805799 became NNA_Mycorrhiza, SRR11805801 became EA_flower, SRR11805803 became NNA_Stem, SRR11805805 became NNA_Flower, SRR11805807 became EA_Mycorrhiza, and SRR11805809 became EA_Stem

#I downloaded the .gff files and put them into my folder using the SFTP window.	
#I have the .gff files, and the reads. BUT I need to pull out the .fasta files from the .gff files. 
#These are two commands to do just that. 
#sed tell it to extract all lines starting from ##FASTA to the end of the files
#The | tells it to send the output of the sed command to the tail command
#tail -n +2 tells it to skip the ##FASTA line, and start at the second line
#the > tells it what to send the files to- in this case, new output .fasta files. 
sed -n '/##FASTA/,$p' EA.gff | tail -n +2 > EA.genes.fasta
sed -n '/##FASTA/,$p' NNA.gff | tail -n +2 > NNA.genes.fasta

#These files are now our references! Yay!
#Now I need the annotations. Those are in the two .tab files which I also downloaded and pulled into the directory from the SFTP window. 

#grep is used for finding specific keywords ina file and filtering, and extracting. It's real good at what it does.
grep 'GO:0015979' EA_mapman_kaas_trinotate.tsv  | cut -f 1 | tr -d '"' > EAphotosynthesis_gene_ids.txt
grep 'GO:0015979' NNA_mapman_kaas_trinotate.tsv  | cut -f 1 | tr -d '"' > NNAphotosynthesis_gene_ids.txt
grep 'GO:0015995' EA_mapman_kaas_trinotate.tsv  | cut -f 1 | tr -d '"' > EAchlorophyll_gene_ids.txt
grep 'GO:0015995' NNA_mapman_kaas_trinotate.tsv  | cut -f 1 | tr -d '"' > NNAchlorophyll_gene_ids.txt

# Define the array for the first species
EA_files=("EA_Flower" "EA_Mycorrhiza" "EA_Stem")

# Cutadapt to trim with --nextseq-trim=30
for EA_ in "${EA_files[@]}"; do
  echo "Cutadapting $EA_ with $threads..."

    # Construct input and output file names
  input_1="${EA_}_1.fastq"
  input_2="${EA_}_2.fastq"
  output_1="${EA_}_1.trimmed.fastq"
  output_2="${EA_}_2.trimmed.fastq"
	
	   # Run cutadapt
    cutadapt --cores 8 --nextseq-trim=30 --minimum-length 60 \
        -o "$output_1" -p "$output_2" \
        "$input_1" "$input_2"
done

# Define the array for the second species
NNA_files=("NNA_Flower" "NNA_Mycorrhiza" "NNA_Stem")

# Cutadapt to trim with --nextseq-trim=30
for NNA_ in "${NNA_files[@]}"; do
    echo "Cutadapting $NNA_ with $threads..."

    # Construct input and output file names
    input_1="${NNA_}_1.fastq"
    input_2="${NNA_}_2.fastq"
    output_1="${NNA_}_1.trimmed.fastq"
    output_2="${NNA_}_2.trimmed.fastq"

    # Run cutadapt
    cutadapt --cores 8 --nextseq-trim=30 --minimum-length 60 \
        -o "$output_1" -p "$output_2" \
        "$input_1" "$input_2"
done
#FASTQC on all reads just for fun (and to take a look at the quality before I go further)

#fastqc -t 8 EA_Stem_1.trimmed.fastq EA_Stem_2.trimmed.fastq
#fastqc -t 8 EA_Flower_1.trimmed.fastq EA_Flower_2.trimmed.fastq
#fastqc -t 8 EA_Mycorrhiza_1.trimmed.fastq EA_Mycorrhiza_2.trimmed.fastq

#fastqc -t 8 NAA_Stem_1.trimmed.fastq NAA_Stem_2.trimmed.fastq
#fastqc -t 8 NAA_Flower_1.trimmed.fastq NAA_Flower_2.trimmed.fastq
#fastqc -t 8 NAA_Mycorrhiza_1.trimmed.fastq NAA_Mycorrhiza_2.trimmed.fastq

echo "Bing bang boom"

refEA=/home/loeffler/Class/Project/Refs/EA.genes.fasta
refNNA=/home/loeffler/Class/Project/Refs/NNA.genes.fasta

#Set a for loop to bbmap all of the lil guys
for EA_ in "${EA_files[@]}"; do
	echo "bbmapping $EA_ with ${threads} threads at minid 95%"
	bbmap.sh ref=$refEA \
	in=${EA_}_1.trimmed.fastq out=mapped_${EA_}_1.sam \
	minid=95 threads=8 rpkm=${EA_}_1_rpkm.txt
	done
echo "All first half reads mapped"

for EA_ in "${EA_files[@]}"; do
	echo "bbmapping $EA_ with ${threads} threads at minid 95%"
	bbmap.sh ref=$refEA \
	in=${EA_}_2.trimmed.fastq out=mapped_${EA_}_2.sam \
	minid=95 threads=8 rpkm=${EA_}_2_rpkm.txt
	done
echo "All second half reads mapped"
echo "All EA reads mapped"
#Same but for NNA
for NNA_ in "${NNA_files[@]}"; do
	echo "bbmapping $NNA_ with ${threads} threads at minid 95%"
	bbmap.sh ref=$refNNA \
	in=${NNA_}_1.trimmed.fastq out=mapped_${NNA_}_1.sam \
	minid=95 threads=8 rpkm=${NNA_}_1_rpkm.txt
	done
echo "All first half reads mapped"

for NNA_ in "${NNA_files[@]}"; do
	echo "bbmapping $NNA_ with ${threads} threads at minid 95%"
	bbmap.sh ref=$refNNA \
	in=${NNA_}_2.trimmed.fastq out=mapped_${NNA_}_2.sam \
	minid=95 threads=8 rpkm=${NNA_}_2_rpkm.txt
	done
echo "All second half reads mapped"
echo "All NNA reads mapped"

#I have the photosynthesis gene IDs in a .txt file already. 
#EAphotosynthesis_gene_ids.txt
#NNAphotosynthesis_gene_ids.txt

#same with the chlorophyll synthesis genes
#EAchlorophyll_gene_ids.txt
#NNAchlorophyll_gene_ids.txt

#Now I just need to go into each of the rpkm.txt files (that's my output from bbmap) and pull out only the photosynthetic genes while keeping the headers so I can tell R to make pretty heatmaps for me.
#grep -Ff mapman_codes.txt input_file.tsv > output_file.tsv

for EA_ in "${EA_files[@]}"; do
	 { head -n 5 "${EA_}_1_rpkm.txt"; grep -Ff EAphotosynthesis_gene_ids.txt "${EA_}_1_rpkm.txt"; } > "${EA_}_pGOI_1.tsv"
	done
echo "EA1 grepped"

for EA_ in "${EA_files[@]}"; do
	 { head -n 5 "${EA_}_2_rpkm.txt"; grep -Ff EAphotosynthesis_gene_ids.txt "${EA_}_2_rpkm.txt"; } > "${EA_}_pGOI_2.tsv"
	done
echo "EA2 grepped"

for NNA_ in "${NNA_files[@]}"; do
	 { head -n 5 "${NNA_}_1_rpkm.txt"; grep -Ff NNAphotosynthesis_gene_ids.txt "${NNA_}_1_rpkm.txt"; } > "${NNA_}_pGOI_1.tsv"
	done
echo "NNA1 grepped"

for NNA_ in "${NNA_files[@]}"; do
	 { head -n 5 "${NNA_}_2_rpkm.txt"; grep -Ff NNAphotosynthesis_gene_ids.txt "${NNA_}_2_rpkm.txt"; } > "${NNA_}_pGOI_2.tsv"
	done
echo "NNA2 grepped"
##and now for the chlorophyll producing genes
for EA_ in "${EA_files[@]}"; do
	 { head -n 5 "${EA_}_1_rpkm.txt"; grep -Ff EAchlorophyll_gene_ids.txt "${EA_}_1_rpkm.txt"; } > "${EA_}_cGOI_1.tsv"
	done
echo "EA1 grepped"

for EA_ in "${EA_files[@]}"; do
	 { head -n 5 "${EA_}_2_rpkm.txt"; grep -Ff EAchlorophyll_gene_ids.txt "${EA_}_2_rpkm.txt"; } > "${EA_}_cGOI_2.tsv"
	done
echo "EA2 grepped"

for NNA_ in "${NNA_files[@]}"; do
	 { head -n 5 "${NNA_}_1_rpkm.txt"; grep -Ff NNAchlorophyll_gene_ids.txt "${NNA_}_1_rpkm.txt"; } > "${NNA_}_cGOI_1.tsv"
	done
echo "NNA1 grepped"

for NNA_ in "${NNA_files[@]}"; do
	 { head -n 5 "${NNA_}_2_rpkm.txt"; grep -Ff NNAchlorophyll_gene_ids.txt "${NNA_}_2_rpkm.txt"; } > "${NNA_}_cGOI_2.tsv"
	done
echo "NNA2 grepped"
##combine the genes of interest files

# Extract relevant columns from both files using a for loop
# Sort and extract RPKM and FPKM values
for EA_ in "${EA_files[@]}"; do
	cut -f1,6,8 "${EA_}_pGOI_1.tsv" | sort -k1,1 > "${EA_}_p_rpkm_fpkm_sorted1.tsv"
	cut -f1,6,8 "${EA_}_pGOI_2.tsv" | sort -k1,1 > "${EA_}_p_rpkm_fpkm_sorted2.tsv"
done
echo "EA files sorted"

for NNA_ in "${NNA_files[@]}"; do
	cut -f1,6,8 "${NNA_}_pGOI_1.tsv" | sort -k1,1 > "${NNA_}_p_rpkm_fpkm_sorted1.tsv"
	cut -f1,6,8 "${NNA_}_pGOI_2.tsv" | sort -k1,1 > "${NNA_}_p_rpkm_fpkm_sorted2.tsv"
done
echo "NNA files sorted"
##Same but for chlorophyll synthesis genes
# Sort and extract RPKM and FPKM values
for EA_ in "${EA_files[@]}"; do
	cut -f1,6,8 "${EA_}_cGOI_1.tsv" | sort -k1,1 > "${EA_}_c_rpkm_fpkm_sorted1.tsv"
	cut -f1,6,8 "${EA_}_cGOI_2.tsv" | sort -k1,1 > "${EA_}_c_rpkm_fpkm_sorted2.tsv"
done
echo "EA files sorted"

for NNA_ in "${NNA_files[@]}"; do
	cut -f1,6,8 "${NNA_}_cGOI_1.tsv" | sort -k1,1 > "${NNA_}_c_rpkm_fpkm_sorted1.tsv"
	cut -f1,6,8 "${NNA_}_cGOI_2.tsv" | sort -k1,1 > "${NNA_}_c_rpkm_fpkm_sorted2.tsv"
done
echo "NNA files sorted"

# Combine the sorted files
for EA_ in "${EA_files[@]}"; do
	join -t $'\t' -1 1 -2 1 "${EA_}_p_rpkm_fpkm_sorted1.tsv" "${EA_}_p_rpkm_fpkm_sorted2.tsv" > "${EA_}_p_combined_temp.tsv"
done

for NNA_ in "${NNA_files[@]}"; do
	join -t $'\t' -1 1 -2 1 "${NNA_}_p_rpkm_fpkm_sorted1.tsv" "${NNA_}_p_rpkm_fpkm_sorted2.tsv" > "${NNA_}_p_combined_temp.tsv"
done
echo "All files combined"
#Do the same with the chlorophyll genes
for EA_ in "${EA_files[@]}"; do
	join -t $'\t' -1 1 -2 1 "${EA_}_c_rpkm_fpkm_sorted1.tsv" "${EA_}_c_rpkm_fpkm_sorted2.tsv" > "${EA_}_c_combined_temp.tsv"
done

for NNA_ in "${NNA_files[@]}"; do
	join -t $'\t' -1 1 -2 1 "${NNA_}_c_rpkm_fpkm_sorted1.tsv" "${NNA_}_c_rpkm_fpkm_sorted2.tsv" > "${NNA_}_c_combined_temp.tsv"
done
echo "All files combined"
# Add headers to the output files
for EA_ in "${EA_files[@]}"; do
	echo -e "#Name\tRPKM_Read1\tFPKM_Read1\tRPKM_Read2\tFPKM_Read2" > "${EA_}_p_combined.tsv"
	cat "${EA_}_p_combined_temp.tsv" >> "${EA_}_p_combined.tsv"
done

for NNA_ in "${NNA_files[@]}"; do
	echo -e "#Name\tRPKM_Read1\tFPKM_Read1\tRPKM_Read2\tFPKM_Read2" > "${NNA_}_p_combined.tsv"
	cat "${NNA_}_p_combined_temp.tsv" >> "${NNA_}_p_combined.tsv"
done
echo "Headers added"
#clean up intermediate files
for EA_ in "${EA_files[@]}"; do
	rm "${EA_}_rpkm_fpkm_sorted1.tsv" "${EA_}_rpkm_fpkm_sorted2.tsv" "${EA_}_p_combined_temp.tsv"
done

for NNA_ in "${NNA_files[@]}"; do
	rm "${NNA_}_rpkm_fpkm_sorted1.tsv" "${NNA_}_rpkm_fpkm_sorted2.tsv" "${NNA_}_p_combined_temp.tsv"
done
echo "Intermediate files cleaned up"
####same with chlorophyll
for EA_ in "${EA_files[@]}"; do
	echo -e "#Name\tRPKM_Read1\tFPKM_Read1\tRPKM_Read2\tFPKM_Read2" > "${EA_}_c_combined.tsv"
	cat "${EA_}_c_combined_temp.tsv" >> "${EA_}_c_combined.tsv"
done

for NNA_ in "${NNA_files[@]}"; do
	echo -e "#Name\tRPKM_Read1\tFPKM_Read1\tRPKM_Read2\tFPKM_Read2" > "${NNA_}_c_combined.tsv"
	cat "${NNA_}_c_combined_temp.tsv" >> "${NNA_}_c_combined.tsv"
done
echo "Headers added"
# Cleanup intermediate files
for EA_ in "${EA_files[@]}"; do
	rm "${EA_}_rpkm_fpkm_sorted1.tsv" "${EA_}_rpkm_fpkm_sorted2.tsv" "${EA_}_c_combined_temp.tsv"
done

for NNA_ in "${NNA_files[@]}"; do
	rm "${NNA_}_rpkm_fpkm_sorted1.tsv" "${NNA_}_rpkm_fpkm_sorted2.tsv" "${NNA_}_c_combined_temp.tsv"
done
echo "Intermediate files cleaned up"


#Now that I have these files, I can go into R to create my figures. I will also provide my .rmd file. 
#here's a stupid thing I did. I was dumb and kept the headers in the combined.tsv files. So the problem is that R won't like to read those, especially since they have a # in front of them. 
#Again, could I have written something to fix this? probably. But I didn't have that many files so I just brute forced it and deleted all but one of the headers and the # in front of 'name' manually. 
#The header to keep is the one that gives titles to all of the columns. That's the only one that matters. 


