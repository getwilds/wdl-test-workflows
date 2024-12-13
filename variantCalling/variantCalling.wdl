version 1.0
## Consensus variant calling workflow for human panel-based DNA sequencing.
## Input requirements:
## - Pair-end sequencing data in unmapped BAM (uBAM) format that comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output :
## - recalibrated bam and it's index and md5
## - GATK vcf
## - Annovar annotated vcfs and tabular file
## 

#### STRUCT DEFINITIONS

# struct for the input files for a given sample
struct SampleInputs {
  String sample_name # Name of the sample in question
  File bam_file # Bam alignment file for the sample in question
  File bed_file # Bed file containing the intervals over which to call variants
}

# struct for all the reference data needed for the run
struct ReferenceData {
  String ref_name # Build version being used as the reference genome, e.g. "hg19", "hg38", etc.
  File ref_fasta # Reference genome fasta file
  File ref_fasta_index # Index for reference genome fasta file
  File ref_dict # Reference genome dictionary file
  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy
  # references such as b37 and hg19.
  File? ref_alt 
  File ref_amb # Reference .amb file for bwa
  File ref_ann # Reference .ann file for bwa
  File ref_bwt # Reference .bwt file for bwa
  File ref_pac # Reference .pac file for bwa
  File ref_sa # Reference .sa file for bwa
  File dbSNP_vcf # dbSNP vcf file for reference during variant calling
  File dbSNP_vcf_index # Index for the dbSNP vcf file
  Array[File] known_indels_sites_VCFs # vcf's of known polymorphic sites to exclude
  Array[File] known_indels_sites_indices # Indexes for the known polymorphic sites vcf's
  String annovar_protocols # Annotation protocols to use, e.g. refGene, cosmic70, etc.
  String annovar_operation # Operation level for annovar, e.g. "g" = gene, "f" = filter, etc.
}

#### WORKFLOW DEFINITION

workflow PanelBwaGatk4Annovar {
  input {
    Array[SampleInputs] sample_batch
    ReferenceData reference_genome
  }

  String gatk_docker = "getwilds/gatk:4.3.0.0"
  String annovar_docker = "getwilds/annovar:${reference_genome.ref_name}"
  String bwa_docker = "getwilds/bwa:0.7.17"

  scatter (sample in sample_batch){
    File bam = sample.bam_file
    File bed = sample.bed_file
    # Get the basename, i.e. strip the filepath and the extension
    String base_file_name = sample.sample_name + "." + reference_genome.ref_name

    # Prepare bed file and check sorting
    call SortBed {
      input:
        unsorted_bed = bed,
        ref_dict = reference_genome.ref_dict,
        task_docker = gatk_docker
    }

    # convert unmapped bam to fastq
    call SamToFastq {
      input:
        input_bam = bam,
        base_file_name = base_file_name,
        task_docker = gatk_docker
    }

    #  Map reads to reference
    call BwaMem {
      input:
        input_fastq = SamToFastq.output_fastq,
        base_file_name = base_file_name,
        ref_fasta = reference_genome.ref_fasta,
        ref_fasta_index = reference_genome.ref_fasta_index,
        ref_dict = reference_genome.ref_dict,
        ref_alt = reference_genome.ref_alt,
        ref_amb = reference_genome.ref_amb,
        ref_ann = reference_genome.ref_ann,
        ref_bwt = reference_genome.ref_bwt,
        ref_pac = reference_genome.ref_pac,
        ref_sa = reference_genome.ref_sa,
        cpu_needed = 4,
        task_docker = bwa_docker
    }

    # Merge original uBAM and BWA-aligned BAM
    call MergeBamAlignment {
      input:
        unmapped_bam = bam,
        aligned_bam = BwaMem.output_bam,
        base_file_name = base_file_name,
        ref_fasta = reference_genome.ref_fasta,
        ref_fasta_index = reference_genome.ref_fasta_index,
        ref_dict = reference_genome.ref_dict,
        task_docker = gatk_docker
    }

    # Generate the recalibration model by interval
    call ApplyBaseRecalibrator {
      input:
        input_bam = MergeBamAlignment.output_bam,
        input_bam_index = MergeBamAlignment.output_bai,
        base_file_name = base_file_name,
        intervals = SortBed.intervals,
        dbSNP_vcf = reference_genome.dbSNP_vcf,
        dbSNP_vcf_index = reference_genome.dbSNP_vcf_index,
        known_indels_sites_VCFs = reference_genome.known_indels_sites_VCFs,
        known_indels_sites_indices = reference_genome.known_indels_sites_indices,
        ref_dict = reference_genome.ref_dict,
        ref_fasta = reference_genome.ref_fasta,
        ref_fasta_index = reference_genome.ref_fasta_index,
        task_docker = gatk_docker
    }

    # Generate haplotype caller vcf
    call HaplotypeCaller {
      input:
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        input_bam_index = ApplyBaseRecalibrator.recalibrated_bai,
        intervals = SortBed.intervals,
        base_file_name = base_file_name,
        ref_dict = reference_genome.ref_dict,
        ref_fasta = reference_genome.ref_fasta,
        ref_fasta_index = reference_genome.ref_fasta_index,
        dbSNP_vcf = reference_genome.dbSNP_vcf,
        task_docker = gatk_docker
    }

    # Annotate variants
    call annovar {
      input:
        input_vcf = HaplotypeCaller.output_vcf,
        ref_name = reference_genome.ref_name,
        base_file_name = base_file_name,
        annovar_operation = reference_genome.annovar_operation,
        annovar_protocols = reference_genome.annovar_protocols,
        task_docker = annovar_docker
    }
  } # End scatter

  # Outputs that will be retained when execution is complete
  output {
    Array[File] analysis_ready_bam = ApplyBaseRecalibrator.recalibrated_bam 
    Array[File] analysis_ready_bai = ApplyBaseRecalibrator.recalibrated_bai
    Array[File] gatk_vcf = HaplotypeCaller.output_vcf
    Array[File] annotated_vcf = annovar.output_annotated_vcf
    Array[File] annotated_table = annovar.output_annotated_table
  }

  parameter_meta {
    sample_batch: "array of structs pointing to the sample data of interest"
    reference_genome: "struct containing the necessary reference files for alignment"

    analysis_ready_bam: "recalibrated bam file produced during alignment steps"
    analysis_ready_bai: "index file for the recalibrated bam"
    gatk_vcf: "raw vcf that comes straight out of HaplotypeCaller"
    annotated_vcf: "annotated vcf produced by Annovar"
    annotated_table: "tsv file containing the same annotation data"
  }
} # End workflow

#### TASK DEFINITIONS

# Prepare bed file and check sorting
task SortBed {
  input {
    File unsorted_bed
    File ref_dict
    String task_docker
  }

  command <<<
    set -eo pipefail
    echo "Sort bed file"
    sort -k1,1V -k2,2n -k3,3n "~{unsorted_bed}" > sorted.bed
    echo "Transform bed file to intervals list with Picard----------------------------------------"
    gatk --java-options "-Xms4g" \
      BedToIntervalList \
      -I sorted.bed \
      -O sorted.interval_list \
      -SD "~{ref_dict}"
  >>>

  output {
    File intervals = "sorted.interval_list"
    File sorted_bed = "sorted.bed"
  }

  runtime {
    docker: task_docker
  }

  parameter_meta {
    unsorted_bed: "unsorted bed file containing the intervals over which to call variants"
    ref_dict: "reference genome dictionary file"
    task_docker: "Docker container to use during execution"

    intervals: "sorted interval list containing the intervals for variant-calling"
    sorted_bed: "sorted bed file containing the intervals for variant-calling"
  }
}

# Read unmapped BAM, convert to FASTQ
task SamToFastq {
  input {
    File input_bam
    String base_file_name
    String task_docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      SamToFastq \
      --INPUT "~{input_bam}" \
      --FASTQ "~{base_file_name}.fastq" \
      --INTERLEAVE true \
      --INCLUDE_NON_PF_READS true 
  >>>

  output {
    File output_fastq = "~{base_file_name}.fastq"
  }

  runtime {
    docker: task_docker
  }

  parameter_meta {
    input_bam: "unmapped bam containing raw reads"
    base_file_name: "base file name to use when saving the fastq"
    task_docker: "Docker container to use during execution"

    output_fastq: "final converted fastq"
  }
}

# align to genome
task BwaMem {
  input {
    File input_fastq
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    Int cpu_needed
    String task_docker
  }

  command <<<
    set -eo pipefail
    bwa mem \
      -p -v 3 -t ~{cpu_needed} -M \
      "~{ref_fasta}" "~{input_fastq}" > "~{base_file_name}.sam"
    samtools view -1bS -@ ~{cpu_needed - 1} \
      -o "~{base_file_name}.aligned.bam" "~{base_file_name}.sam"
  >>>

  output {
    File output_bam = "~{base_file_name}.aligned.bam"
  }

  runtime {
    docker: task_docker
    memory: "4GB"
    cpu: cpu_needed
  }

  parameter_meta {
    input_fastq: "fastq file containing raw reads to be aligned"
    base_file_name: "base file name to use when saving the bam file"
    ref_fasta: "reference genome fasta file"
    ref_fasta_index: "index for the reference genome fasta file"
    ref_dict: "reference genome dictionary file"
    ref_alt: "reference genome .alt file for bwa (optional)"
    ref_amb: "reference genome .amb file for bwa"
    ref_ann: "reference genome .ann file for bwa"
    ref_bwt: "reference genome .bwt file for bwa"
    ref_pac: "reference genome .pac file for bwa"
    ref_sa: "reference genome .sa file for bwa"
    cpu_needed: "number of cpus to use during alignment"
    task_docker: "Docker container to use during execution"

    output_bam: "bam alignment file containing only aligned reads"
  }
}

# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  input {
    File unmapped_bam
    File aligned_bam
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String task_docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Xmx8g" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM "~{aligned_bam}" \
      --UNMAPPED_BAM "~{unmapped_bam}" \
      --OUTPUT "~{base_file_name}.merged.bam" \
      --REFERENCE_SEQUENCE "~{ref_fasta}" \
      --PAIRED_RUN true \
      --SORT_ORDER coordinate \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 200000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --CREATE_INDEX true
  >>>

  output {
    File output_bam = "~{base_file_name}.merged.bam"
    File output_bai = "~{base_file_name}.merged.bai"
  }

  runtime {
    docker: task_docker
  }

  parameter_meta {
    unmapped_bam: "unmapped bam containing raw reads"
    aligned_bam: "mapped bam containing aligned reads"
    base_file_name: "base file name to use when saving the merged bam file"
    ref_fasta: "reference genome fasta file"
    ref_fasta_index: "index for the reference genome fasta file"
    ref_dict: "reference genome dictionary file"
    task_docker: "Docker container to use during execution"

    output_bam: "merged bam file containing both mapped and unmapped reads"
    output_bai: "index for the merged bam file"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model and apply it
task ApplyBaseRecalibrator {
  input {
    File input_bam
    File intervals 
    File input_bam_index
    String base_file_name
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String task_docker
  }

  command <<<
    set -eo pipefail
    samtools index "~{input_bam}"
    gatk --java-options "-Xms4g" \
        BaseRecalibrator \
        -R "~{ref_fasta}" \
        -I "~{input_bam}" \
        -O "~{base_file_name}.recal_data.csv" \
        --known-sites "~{dbSNP_vcf}" \
        --known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
        --intervals "~{intervals}" \
        --interval-padding 100 
    gatk --java-options "-Xms4g" \
        ApplyBQSR \
        -bqsr "~{base_file_name}.recal_data.csv" \
        -I "~{input_bam}" \
        -O "~{base_file_name}.recal.bam" \
        -R "~{ref_fasta}" \
        --intervals "~{intervals}" \
        --interval-padding 100 
    samtools view -H "~{base_file_name}.recal.bam" | grep @SQ \
      | sed 's/@SQ\tSN:\|LN://g' > "~{base_file_name}.sortOrder.txt"
  >>>

  output {
    File recalibrated_bam = "~{base_file_name}.recal.bam"
    File recalibrated_bai = "~{base_file_name}.recal.bai"
    File sortOrder = "~{base_file_name}.sortOrder.txt"
  }

  runtime {
    docker: task_docker
  }

  parameter_meta {
    input_bam: "merged bam file to undergo recalibration"
    intervals: "interval list containing the regions for recalibration"
    input_bam_index: "index for the merged input bam file"
    base_file_name: "base file name to use when saving the recalibrated bam file"
    dbSNP_vcf: "dbSNP vcf file for reference during recalibration"
    dbSNP_vcf_index: "index for the dbSNP vcf file"
    known_indels_sites_VCFs: "vcf's of known polymorphic sites to exclude"
    known_indels_sites_indices: "indexes for the known polymorphic sites vcf's"
    ref_fasta: "reference genome fasta file"
    ref_fasta_index: "index for the reference genome fasta file"
    ref_dict: "reference genome dictionary file"
    task_docker: "Docker container to use during execution"

    recalibrated_bam: "recalibrated bam file produced by BQSR"
    recalibrated_bai: "index for the recalibrated bam file"
    sortOrder: "text file containing the sorting order if necessary"
  }
}

# HaplotypeCaller per-sample
task HaplotypeCaller {
  input {
    File input_bam
    File input_bam_index
    String base_file_name
    File intervals
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File dbSNP_vcf
    String task_docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Xmx4g" \
      HaplotypeCaller \
      -R "~{ref_fasta}" \
      -I "~{input_bam}" \
      -O "~{base_file_name}.GATK.vcf" \
      --intervals "~{intervals}" \
      --interval-padding 100 
  >>>

  output {
    File output_vcf = "~{base_file_name}.GATK.vcf"
    File output_vcf_index = "~{base_file_name}.GATK.vcf.idx"
  }

  runtime {
    docker: task_docker
  }

  parameter_meta {
    input_bam: "recalibrated bam file to use for variant calling"
    input_bam_index: "index for the recalibrated bam file"
    base_file_name: "base file name to use when saving the vcf"
    intervals: "sorted interval list containing the intervals for variant-calling"
    ref_dict: "reference genome dictionary file"
    ref_fasta: "reference genome fasta file"
    ref_fasta_index: "index for the reference genome fasta file"
    dbSNP_vcf: "dbSNP vcf file for reference during variant calling"
    task_docker: "Docker container to use during execution"

    output_vcf: "resulting vcf containing the raw variant calls from HaplotypeCaller"
    output_vcf_index: "index for the resulting vcf"
  }
}

# annotate with annovar
task annovar {
  input {
    File input_vcf
    String base_file_name
    String ref_name
    String annovar_protocols
    String annovar_operation
    String task_docker
    String base_vcf_name = basename(input_vcf, ".vcf")
  }
  
  command <<<
    set -eo pipefail
    perl /annovar/table_annovar.pl "~{input_vcf}" /annovar/humandb/ \
      -buildver "~{ref_name}" \
      -outfile "~{base_vcf_name}" \
      -remove \
      -protocol "~{annovar_protocols}" \
      -operation "~{annovar_operation}" \
      -nastring . -vcfinput
  >>>

  output {
    File output_annotated_vcf = "~{base_file_name}.GATK.~{ref_name}_multianno.vcf"
    File output_annotated_table = "~{base_file_name}.GATK.~{ref_name}_multianno.txt"
  }

  runtime {
    docker: task_docker
  }

  parameter_meta {
    input_vcf: "vcf containing variant calls to be annotated by Annovar"
    base_file_name: "base file name to use when saving the annotation results"
    ref_name: "Build version being used as the reference genome, e.g. 'hg19', 'hg38', etc."
    annovar_protocols: "Annotation protocols to use, e.g. refGene, cosmic70, etc."
    annovar_operation: "Operation level for annovar, e.g. 'g' = gene, 'f' = filter, etc."
    task_docker: "Docker container to use during execution"
    base_vcf_name: "base vcf name to be used within Annovar"

    output_annotated_vcf: "annotated vcf produced by Annovar"
    output_annotated_table: "tsv file containing the same annotation data"
  }
}