require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/workflow'

module PCAWG
  extend Resource
  self.subdir = 'share/data/projects/PCAWG'

  WGS_SAMPLE_FIELD = 'tumor_wgs_aliquot_id'
  DONOR_FIELD = 'icgc_donor_id'
  RNA_TUMOR_SAMPLE = 'RNA_tumor_sample'
  RNA_NORMAL_SAMPLE = 'RNA_normal_sample'


  DATA_DIR = Rbbt.root.share.data.projects.PCAWG[".source"]

  PROJECT_VAR_DIR = Rbbt.root.var.PCAWG

  def self.organism
    "Hsa/feb2014"
  end

  PCAWG.claim PCAWG.selected_donor_samples, :proc do |filename|
    list = TSV.traverse DATA_DIR["aliquot_donor_tumor.whitelist.tsv.gz"].find, :into => [], :type => :array do |line|
      next if line =~ /DONOR_UNIQ/
      line.split("\t")[1]
    end
    list * "\n" + "\n"
  end

  PCAWG.claim PCAWG.donor_wgs_samples, :proc do |filename|
    tsv = DATA_DIR['release_may2016.v1.2.tsv'].tsv :key_field => WGS_SAMPLE_FIELD, :fields => [DONOR_FIELD], :header_hash => '', :sep2 => ','
    selected_donor_samples = PCAWG.selected_donor_samples.list
    tsv = tsv.select selected_donor_samples
    tsv = tsv.reorder 'icgc_donor_id', %w(tumor_wgs_aliquot_id)
    tsv.to_single.to_s
  end
  
  PCAWG.claim PCAWG.donor_rna_samples, :proc do |filename|
    tsv = DATA_DIR['release_may2016.v1.2.tsv'].tsv :key_field => DONOR_FIELD, :fields => %w(tumor_rna_seq_star_alignment_bam_file_name normal_rna_seq_star_alignment_bam_file_name), :header_hash => '', :sep2 => ','
    tsv.process "tumor_rna_seq_star_alignment_bam_file_name" do |names|
      names.collect{|name|
        name.split(".")[1]
      }
    end
    
    tsv.process "normal_rna_seq_star_alignment_bam_file_name" do |names|
      names.collect{|name|
        name.split(".")[1]
      }
    end

    tsv.fields = [RNA_TUMOR_SAMPLE, RNA_NORMAL_SAMPLE]

    tsv.to_list.to_s
  end

  PCAWG.claim PCAWG.donor_histology, :proc do |filename|
    fields = %w(organ_system histology_abbreviation histology_tier1 histology_tier2 histology_tier3 histology_tier4 tumour_histological_code tumour_histological_type tumour_stage tumour_grade specimen_donor_treatment_type)
    tsv = DATA_DIR['pcawg_specimen_histology_August2016_v1.tsv'].tsv :key_field => 'icgc_donor_id', :fields => fields, :type => :double, :sep2 => /,\s*/
  end

  PCAWG.claim PCAWG.donor_clinical, :proc do |filename|
    fields = %w(donor_sex donor_vital_status donor_diagnosis_icd10 first_therapy_type first_therapy_response donor_age_at_diagnosis donor_survival_time donor_interval_of_last_followup tobacco_smoking_history_indicator tobacco_smoking_intensity alcohol_history alcohol_history_intensity)
    tsv = DATA_DIR["pcawg_donor_clinical_August2016_v1.tsv - pcawg_donor_clinical_August2016_v1.tsv.tsv"].tsv :key_field => 'icgc_donor_id', :fields => fields, :type => :list
    tsv.to_s
  end

  PCAWG.claim DATA_DIR['joint_fpkm_uq.tsv.gz'].tap{|o| o.resource = PCAWG}, :proc do |filename|
    raise "Please place the file joint_fpkm_uq.tsv.gz into #{ filename }"
  end

  PCAWG.claim PCAWG.matrices.gene_expression, :proc do 
    TSV.traverse DATA_DIR["joint_fpkm_uq.tsv.gz"].find, :type => :array, :into => :stream, :bar => true do |line|
      if line =~ /^feature/
        fields = line.split("\t")
        fields[0] = "Ensembl Gene ID"
        "#" + fields * "\t"
      else
        values = line.split("\t")
          values[0] = values[0].split(".").first
          values * "\t"
      end
    end
  end

  PCAWG.claim DATA_DIR['final_consensus_12aug_passonly_whitelist_31aug_snv_indel_v3.maf.gz'].tap{|o| o.resource = PCAWG}, :proc do |filename|
    raise "You do not have permission to view Genomic Mutations in this server. Otherwise, please place the file final_consensus_12aug_passonly_whitelist_31aug_snv_indel_v3.maf.gz into #{ filename }"
  end

  PCAWG.claim PCAWG.genotypes, :proc do |real|
    TmpFile.with_file do |directory|
      Path.setup(directory)
      aliquote2donor = PCAWG.donor_wgs_samples.tsv :key_field => 'tumor_wgs_aliquot_id', :type => :single
      if Open.exists? directory
        FileUtils.rm_rf directory.find
      end
      FileUtils.mkdir_p directory 
      last_donor = nil
      io = nil
      TSV.traverse DATA_DIR['final_consensus_12aug_passonly_whitelist_31aug_snv_indel_v3.maf.gz'], :type => :array, :bar => true do |line|
        next if line =~ /^Tumor_Sample_Barcode/
        parts = line.split("\t")
        ali, chr, start, eend, ref, alt, alt2  = parts.values_at 0,3,4,5,10,11

        donor = aliquote2donor[ali]
        next if donor.nil?
        start = start.to_i
        eend = eend.to_i

        _muts = ([alt, alt2].compact.uniq - [ref])
        _mut  = _muts.first

        pos, muts = Misc.correct_vcf_mutation start, ref, _mut
        mutations = muts.collect{|m| [chr, pos, m] * ":"}

        if last_donor != donor
          io.close if io
          if File.exists?(directory[donor])
            io = Open.open(directory[donor], :mode => 'a')
          else
            io = Open.open(directory[donor], :mode => 'w')
          end
          io.puts(mutations * "\n")
        else
          io.puts(mutations * "\n")
        end
      end
      io.close
      FileUtils.mv directory.find, real.find
    end
    nil
  end
end
