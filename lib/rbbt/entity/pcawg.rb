
require 'rbbt/entity/pcawg/study'
require 'rbbt/entity/pcawg/sample'
require 'rbbt/entity/pcawg/specimen'
require 'rbbt/entity/pcawg/donor'

module PCAWG

  def self.all_donors
    @@all_donors ||= PCAWG.donor_samples.tsv(:key_field => PCAWG::DONOR_FIELD, :fields => [], :type => :list).keys
  end

  def self.sample2donor
    @@sample2donor ||= PCAWG.donor_samples.index :target => PCAWG::DONOR_FIELD, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache
  end

  def self.specimen_histology_index
    @@specimen_histology_index ||= PCAWG.specimen_histology.tsv :key_field => PCAWG::SPECIMEN_FIELD, :fields => ["organ_system", "histology_abbreviation", "histology_tier1", "histology_tier2", "histology_tier3", "histology_tier4", "tumour_stage", "tumour_grade", "percentage_cellularity", "level_of_cellularity", "tcga_expert_re-review", "specimen_donor_treatment_type"], :type => :list
  end
  
  def self.specimen2donor
    @@specimen2donor ||= self.sample2donor
  end

  def self.sample2specimen
    @@sample2specimen ||= PCAWG.specimen_histology.index :fields => ['tcga_sample_uuid'], :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache
  end

  def self.donor2tumor_specimen
    @@donor2tumor_specimen ||= PCAWG.donor_sample_info.tsv :key_field => PCAWG::DONOR_FIELD, :fields => ['tumor_wgs_icgc_specimen_id'], :type => :flat, :merge => true, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache
  end

  def self.donor2normal_specimen
    @@donor2normal_specimen ||= PCAWG.donor_samples.tsv :key_field => PCAWG::DONOR_FIELD, :fields => ['normal_wgs_icgc_specimen_id'], :type => :flat, :merge => true, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache
  end

  def self.donor2SNV_sample
    @@donor2SNV_sample ||= PCAWG.donor_samples.tsv :key_field => PCAWG::DONOR_FIELD, :fields => [PCAWG::SNV_SAMPLE_FIELD], :type => :flat, :merge => true, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache
  end

  def self.donor2CNV_sample
    @@donor2CNV_sample ||= Hash.new{[]}
  end

  def self.donor2expression_samples
    @@donor2expression_samples ||= PCAWG.donor_samples.tsv :key_field => PCAWG::DONOR_FIELD, :fields => ['RNA Seq ID'], :type => :flat, :merge => true, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache
  end

end
