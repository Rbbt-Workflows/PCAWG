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

  def self.abbr_color(abbr)
    @@abbr_color ||= PCAWG.subtype_info.tsv :fields => ["Color (RGB code)"], :type => :single
    @@abbr_color[abbr]
  end

  PCAWG.claim PCAWG.subtype_info, :proc do
    file = PCAWG::DATA_DIR["tumour_subtype_consolidation_map.tsv - Unique List of Tumour Types_August.tsv"]
    tsv = file.tsv :type => :list, :key_field => "Abbreviation", :header_hash => '', :grep => "^[[:space:]]\|MISSING", :invert_grep => true, :fix => Proc.new{|l| l.gsub(/([a-z]+)CA(\s)/, '\1Ca\2')}
    tsv.to_s
  end

  PCAWG.claim PCAWG.enhancer_ranges, :proc do
    file = PCAWG::DATA_DIR.annotations["map.enhancer.gene"]
    TSV.traverse file, :type => :array, :into => :stream do |line|
      range, genes = line.split("\t")
      chr,start,eend = range.split(/:|-/)
      chr.sub!('chr','')
      [chr, start, eend, genes.split(";").uniq*";"] * "\t"
    end
  end

  PCAWG.claim PCAWG.selected_donor_samples, :proc do |filename|
    list = TSV.traverse DATA_DIR["aliquot_donor_tumor.whitelist.tsv.gz"].find, :into => [], :type => :array do |line|
      next if line =~ /DONOR_UNIQ/
      line.split("\t")[1]
    end
    list * "\n" + "\n"
  end

  PCAWG.claim PCAWG.donor_specimens, :proc do |filename|
    tsv = DATA_DIR['release_may2016.v1.2.tsv'].tsv :key_field => DONOR_FIELD, :fields => %w(normal_wgs_icgc_specimen_id tumor_wgs_icgc_specimen_id normal_rna_seq_icgc_specimen_id tumor_rna_seq_icgc_specimen_id), :header_hash => '', :sep2 => ','
    tsv.fields = ["Normal WGS Specimen", "Tumor WGS Specimen", "Normal RNA Specimen", "Tumor RNA Specimen"]
    tsv.to_s
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
    tsv = DATA_DIR['pcawg_specimen_histology_August2016_v3.tsv'].tsv :key_field => 'icgc_donor_id', :fields => fields, :type => :double, :sep2 => /,\s*/, :fix => Proc.new {|l| l.gsub(/([a-z]+)CA(\s)/, '\1Ca\2')}
  end

  PCAWG.claim PCAWG.donor_clinical, :proc do |filename|
    fields = %w(donor_sex donor_vital_status donor_diagnosis_icd10 first_therapy_type first_therapy_response donor_age_at_diagnosis donor_survival_time donor_interval_of_last_followup tobacco_smoking_history_indicator tobacco_smoking_intensity alcohol_history alcohol_history_intensity)
    tsv = DATA_DIR["pcawg_donor_clinical_August2016_v3.tsv"].tsv :key_field => 'icgc_donor_id', :fields => fields, :type => :list, :fix => Proc.new{|l| l.gsub(/([a-z]+)CA(\s)/, '\1Ca\2')}
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

  PCAWG.claim PCAWG.genotypes_random_Inigo, :proc do |real|
    TmpFile.with_file do |directory|
      Path.setup(directory)
      aliquote2donor = PCAWG.donor_wgs_samples.tsv :key_field => 'tumor_wgs_aliquot_id', :type => :single
      if Open.exists? directory
        FileUtils.rm_rf directory.find
      end
      FileUtils.mkdir_p directory 
      last_donor = nil
      io = nil
      TSV.traverse DATA_DIR['Inigo_randomized_neutral.6cols'], :type => :array, :bar => true do |line|
        next if line =~ /^Tumor_Sample_Barcode/
        parts = line.split("\t")
        ali, chr, start, ref, alt  = parts
        alt2 = nil

        donor = aliquote2donor[ali]
        next if donor.nil?
        start = start.to_i

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

  PCAWG.claim DATA_DIR['AugustRelease_Simulations_Broad.maf'].tap{|o| o.resource = PCAWG}, :proc do |filename|
    raise "Please place the file AugustRelease_Simulations_Broad.maf.gz into #{ filename }"
  end

  PCAWG.claim PCAWG.genotypes_random_Broad, :proc do |real|
    TmpFile.with_file do |directory|
      Path.setup(directory)
      aliquote2donor = PCAWG.donor_wgs_samples.tsv :key_field => 'tumor_wgs_aliquot_id', :type => :single
      if Open.exists? directory
        FileUtils.rm_rf directory.find
      end
      FileUtils.mkdir_p directory 
      last_donor = nil
      io = nil
      TSV.traverse DATA_DIR['AugustRelease_Simulations_Broad.maf'], :type => :array, :bar => true do |line|
        next if line =~ /^Tumor_Sample_Barcode/
        parts = line.split("\t")
        chr, start, ali, ref, alt  = parts
        alt2 = nil

        donor = aliquote2donor[ali]
        next if donor.nil?
        start = start.to_i

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

  PCAWG.claim DATA_DIR['final_consensus_12aug_passonly_whitelist_31aug_snv_indel_v3.dkfz_randomize_25kbwindow_complete_set.maf'].tap{|o| o.resource = PCAWG}, :proc do |filename|
    raise "Please place the file final_consensus_12aug_passonly_whitelist_31aug_snv_indel_v3.dkfz_randomize_25kbwindow_complete_set.maf into #{ filename }"
  end

  PCAWG.claim PCAWG.genotypes_random_DKFZ, :proc do |real|
    TmpFile.with_file do |directory|
      Path.setup(directory)
      aliquote2donor = PCAWG.donor_wgs_samples.tsv :key_field => 'tumor_wgs_aliquot_id', :type => :single
      if Open.exists? directory
        FileUtils.rm_rf directory.find
      end
      FileUtils.mkdir_p directory 
      last_donor = nil
      io = nil
      TSV.traverse DATA_DIR['final_consensus_12aug_passonly_whitelist_31aug_snv_indel_v3.dkfz_randomize_25kbwindow_complete_set.maf'], :type => :array, :bar => true do |line|
        next if line =~ /^Tumor_Sample_Barcode/
        parts = line.split("\t")
        ali, chr, start, eend, ref, alt, alt2  = parts.values_at 0, 3,4,5,9,10,11

        donor = aliquote2donor[ali]
        next if donor.nil?
        start = start.to_i

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

  PCAWG.claim DATA_DIR['October_2016_whitelist_2583.snv_mnv_indel.maf.gz'].tap{|o| o.resource = PCAWG}, :proc do |filename|
    raise "You do not have permission to view Genomic Mutations in this server. Otherwise, please place the file October_2016_whitelist_2583.snv_mnv_indel.maf.gz into #{ filename }"
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
      TSV.traverse DATA_DIR['October_2016_whitelist_2583.snv_mnv_indel.maf.gz'], :type => :array, :bar => true do |line|
        next if line =~ /^Tumor_Sample_Barcode/
        parts = line.split("\t")
        chr, start, eend, ref, alt, alt2, ali  = parts.values_at 1,2,3,7,8,9,12

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

  PCAWG.claim PCAWG.SV, :proc do |real|
    file = DATA_DIR['pcawg_consensus_1.5.160912.somatic.sv.tar.gz'].find
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        directory.glob("*.bedpe").each do |file|
          sample = File.basename(file).split(".").first
          donor = sample2donor[sample]
          next if donor.nil?
          stream = TSV.traverse Open.open(file), :type => :array, :into => :stream do |line|
            next if line.include? 'chrom'
            chr1, pos1, _pos1, chr2, pos2, _pos2, _id, support, strand1, strand2, sv_class = line.split("\t")
            [chr1, pos1, sv_class,chr2,pos2] * ":"
          end

          Open.write(real[donor].find, stream.read)
        end
      end
    end
    FileUtils.rm_rf directory
    nil
  end
end
