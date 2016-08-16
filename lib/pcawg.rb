require 'rbbt/resource'
module PCAWG
  extend Resource
  self.subdir = 'share/data/projects/PCAWG'

  PROJECT_DIR= PCAWG.final

  def self.organism
    "Hsa/feb2014"
  end

  PCAWG.claim PROJECT_DIR.genotypes, :proc do |dir|
    tar = PCAWG.preliminary_final_release['preliminary_final_release.snvs.tgz'].find

    FileUtils.mkdir_p dir unless File.directory? dir
    TmpFile.with_file do |tmpdir|
      Path.setup(tmpdir)
      Misc.in_dir tmpdir do
        CMD.cmd("tar xvfz #{tar}")

        tmpdir.preliminary_final_release.snv_mnv.glob('*.tbi').each{|f| FileUtils.rm f}

        tmpdir.preliminary_final_release.snv_mnv.glob('*.gz').each do |f| 
          name = File.basename(f).sub('.annotated.snv_mnv.vcf.gz','.vcf')
          CMD.cmd("zcat #{f} | grep -v LOWSUPPORT > #{dir[name]}")
        end

        dir = Path.setup(dir)
        dir.glob('*.vcf').each do |f|
          CMD.cmd("gzip #{f}")
        end

      end
    end

  end

  PCAWG.claim PCAWG.donor_samples, :proc do
    release = PCAWG.preliminary_final_release["release_may2016.v1.tsv"].tsv :header_hash => "", :type => :double, :key_field => 'icgc_donor_id', :sep2 => ','

    release.add_field "RNA Seq ID" do |key,values|
      values[81].collect{|id| id.split('.')[1]}
    end

    release.to_s
  end

  PCAWG.claim PCAWG.specimen_histology, :proc do
    histology = PCAWG.preliminary_final_release["pcawg_specimen_histology_May2016_v2.tsv"].tsv :type => :double, :key_field => 'icgc_specimen_id', :sep2 => ','
    histology.to_s
  end

  PCAWG.claim PCAWG.subtype_info, :proc do
    tsv = PCAWG.preliminary_final_release["tumour_subtype_consolidation_map.tsv - Unique List of Tumour Types_May.tsv"].tsv :type => :list, :key_field => "Abbreviation"
    ppp tsv.to_s
    tsv.to_s
  end

  def self.all_abbreviations
    PCAWG.specimen_histology.tsv(:key_field => "histology_abbreviation", :fields => []).keys
  end

  def self.all_histologies(tier = 1)
    PCAWG.specimen_histology.tsv(:key_field => "histology_tier"+tier.to_s, :fields => []).keys
  end

  def self.specimens_with_histology(h, key = nil)
    key, h = h.split ":" if h.include?(":") and key.nil?
    key ||= 'histology_abbreviation'
    PCAWG.specimen_histology.tsv.select(key => h).keys
  end

  def self.specimen_histology_key(s, key = 'histology_abbreviation')
    @@speciment_histology_index ||= {}
    @@speciment_histology_index[key] ||= PCAWG.specimen_histology.tsv(:fields => [key], :key_field => "icgc_specimen_id", :type => :single)
    @@speciment_histology_index[key][s]
  end

  def self.specimen_donor(specimen)
    @@index ||= PCAWG.donor_samples.index :target => 'icgc_donor_id'
    @@index[specimen]
  end

  def self.sample_donor(sample)
    @@index ||= PCAWG.donor_samples.index :target => 'icgc_donor_id'
    @@index[sample]
  end

  def self.donors_with_histology(*args)
    specimens = specimens_with_histology(*args)
    donors = specimens.collect{|specimen| PCAWG.specimen_donor(specimen)}.uniq
    Donor.setup(donors)
    donors.extend AnnotatedArray
    donors
  end

  def self.donor_histology(d, *args)
    specimens = Donor.donor2tumor_specimen[d]
    return specimens.collect{|s| specimen_histology_key(s, *args)}.compact.first
  end


  def self.abb_color(abb)
    @@color_index ||= PCAWG.subtype_info.tsv :fields => ["Color (RGB code)"], :type => :single, :key_field => "Abbreviation"
    @@color_index[abb]
  end

end

Workflow.require_workflow "Sample"

module Sample
  helper :watson do
    true
  end

  helper :organism do
    PCAWG.organism
  end

  task :organism => :string do
    PCAWG.organism
  end

  task :genomic_mutations => :array do
    raise "Genomic Mutations not accessible for #{clean_name} due to access limitations"
  end
end

module Donor
  extend Entity

  def self.all_donors
    PCAWG.donor_samples.tsv(:key_field => 'icgc_donor_id', :fields => [], :type => :list).keys
  end

  def self.donor2tumor_specimen
    PCAWG.donor_samples.tsv :key_field => 'icgc_donor_id', :fields => ['tumor_wgs_icgc_specimen_id'], :type => :flat, :merge => true, :persist => true, :persist_dir => Rbbt.var.PCAWG.find(:lib)
  end

  def self.donor2normal_specimen
    PCAWG.donor_samples.tsv :key_field => 'icgc_donor_id', :fields => ['normal_wgs_icgc_specimen_id'], :type => :flat, :merge => true, :persist => true, :persist_dir => Rbbt.var.PCAWG.find(:lib)
  end

  def self.donor2SNV_sample
    PCAWG.donor_samples.tsv :key_field => 'icgc_donor_id', :fields => ['dkfz_embl_variant_calling_file_name_prefix'], :type => :flat, :merge => true, :persist => true, :persist_dir => Rbbt.var.PCAWG.find(:lib)
  end

  property :histology_abbreviation => :single do
    specimens = Donor.donor2tumor_specimen[self]
    index = PCAWG.specimen_histology.tsv :key_field => 'icgc_specimen_id', :fields => ['histology_abbreviation'], :type => :flat, :merge => true, :persist => true, :persist_dir => Rbbt.var.PCAWG.find(:lib)
    index.chunked_values_at(specimens).flatten
  end

  property :tumor_specimens => :single do
    Donor.donor2tumor_specimen[self]
  end

  property :normal_specimens => :single do
    Donor.donor2normal_specimen[self]
  end

  property :SNV_sample => :single do
    Sample.setup(Donor.donor2SNV_sample[self], :cohort => 'PCAWG').extend AnnotatedArray
  end

  property :abbr_color => :single do
    abbr = PCAWG.donor_histology(self, 'histology_abbreviation')
    PCAWG.abb_color(abbr)
  end

  property :gene_status => :single2array do |genes|
    gms = self.SNV_sample.first.gene_mutation_status
    gms.chunked_values_at(genes).collect do |values|
      if values
        case 
        when values[:broken] == 'true'
          :broken
        when values[:affected] == 'true'
          :affected
        when values[:affected] == 'loss'
          :loss
        when values[:affected] == 'gain'
          :gain
        end
      else
        nil
      end
    end

  end
end

Workflow.require_workflow "Study"

module Study

  helper :organism do
    PCAWG.organism
  end

  task :organism => :string do
    PCAWG.organism
  end

  task :mappable_genes => :array do
    []
  end

  property :organism => :single do
    PCAWG.organism
  end

  property :donors => :single do
    PCAWG.donors_with_histology(self)
  end

  property :genotyped_samples => :single do
    samples = Sample.setup(donors.collect{|donor| donor.SNV_sample}.compact)
    samples.extend AnnotatedArray
    samples.cohort = "PCAWG"
    samples
  end

  property :has_cnv? => :single do
    false
  end

end

allowed_tasks = %w(mi mi_truncated mi_damaged firestar gene_mutation_status sample_genes)
Study.tasks.keys.each do |task|
  next unless allowed_tasks.include? task.to_s

  Study.export task
end

Sample.tasks.keys.each do |task|
  next unless allowed_tasks.include? task.to_s

  Sample.export task
end
