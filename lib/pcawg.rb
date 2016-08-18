require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/workflow'

module PCAWG
  extend Resource
  self.subdir = 'share/data/projects/PCAWG'

  PROJECT_DIR= PCAWG.final

  DATA_DIR=Rbbt.data.final

  PROJECT_VAR_DIR = Rbbt.var.PCAWG

  SNV_SAMPLE_FIELD_ORIG = 'tumor_wgs_aliquot_id'
  EXPRESSION_SAMPLE_FIELD_ORIG = 'tumor_rna_seq_star_alignment_bam_file_name'
  NORMAL_EXPRESSION_SAMPLE_FIELD_ORIG = 'normal_rna_seq_star_alignment_bam_file_name'

  DONOR_FIELD = 'icgc_donor_id'
  SPECIMEN_FIELD = 'icgc_specimen_id'
  SAMPLE_FIELD = 'icgc_sample_id'

  SNV_SAMPLE_FIELD = 'SNV ID'
  EXPRESSION_SAMPLE_FIELD = 'RNA Seq ID'
  NORMAL_EXPRESSION_SAMPLE_FIELD = 'Normal RNA Seq ID'

  def self.organism
    "Hsa/feb2014"
  end

  #PCAWG.claim PROJECT_DIR.genotypes, :proc do |dir|
  #  tar = PCAWG.preliminary_final_release['final_consensus_12aug_passonly_snv_indel.tgz'].find

  #  FileUtils.mkdir_p dir unless File.directory? dir
  #  TmpFile.with_file do |tmpdir|
  #    Path.setup(tmpdir)
  #    Misc.in_dir tmpdir do
  #      CMD.cmd("tar xvfz #{tar}")

  #      tmpdir.preliminary_final_release.snv_mnv.glob('*.tbi').each{|f| FileUtils.rm f}

  #      tmpdir.preliminary_final_release.snv_mnv.glob('*.gz').each do |f| 
  #        name = File.basename(f).sub('.annotated.snv_mnv.vcf.gz','.vcf')
  #        CMD.cmd("zcat #{f} | grep -v LOWSUPPORT > #{dir[name]}")
  #      end

  #      dir = Path.setup(dir)
  #      dir.glob('*.vcf').each do |f|
  #        CMD.cmd("gzip #{f}")
  #      end

  #    end
  #  end

  #end
  PCAWG.claim PCAWG.blacklisted_donors, :proc do
    tsv = PCAWG::DATA_DIR["release_may2016.v1.blacklisted_donors.tsv"].tsv :header_hash => "", :type => :double, :key_field => DONOR_FIELD, :sep2 => ',', :fields => []
    tsv.keys * "\n"
  end

  PCAWG.claim PCAWG.donor_sample_info, :proc do
    release = PCAWG::DATA_DIR["release_may2016.v1.tsv"].tsv :header_hash => "", :type => :double, :key_field => DONOR_FIELD, :sep2 => ','
    release_BL = PCAWG::DATA_DIR["release_may2016.v1.blacklisted_donors.tsv"].tsv :header_hash => "", :type => :double, :key_field => DONOR_FIELD, :sep2 => ','

    release.merge!(release_BL)

    release = release.add_field EXPRESSION_SAMPLE_FIELD do |key,values|
      values[EXPRESSION_SAMPLE_FIELD_ORIG].collect{|id| id.split('.')[1]}
    end

    release = release.add_field NORMAL_EXPRESSION_SAMPLE_FIELD do |key,values|
      values[NORMAL_EXPRESSION_SAMPLE_FIELD_ORIG].collect{|id| id.split('.')[1]}
    end

    fields = release.fields
    specimen_fields = fields.select{|f| f.include? 'icgc_specimen'}
    specimen_field_positions = specimen_fields.collect{|f| fields.index f}
    release = release.add_field SPECIMEN_FIELD do |key,values|
      values.values_at(*specimen_field_positions).flatten.compact.reject{|s| s.empty? }
    end

    fields = release.fields
    sample_fields = fields.select{|f| f.include? 'icgc_sample'}
    sample_field_positions = sample_fields.collect{|f| fields.index f}
    release = release.add_field SAMPLE_FIELD do |key,values|
      values.values_at(*sample_field_positions).flatten.compact.reject{|s| s.empty? }
    end

    release.to_s
  end
  
  PCAWG.claim PCAWG.preferred_samples, :proc do
    donor = PCAWG.donor_sample_info.tsv :fields => 'donor_unique_id', :key_field => PCAWG::DONOR_FIELD, :type => :list
    tsv = PCAWG::DATA_DIR["PCAWG multi-tumour list - Selection of representative aliquots (Aug 18).tsv"].tsv :header_hash => "", :type => :list, :fields => ["tumor_wgs_aliquot_id"] 
    Log.tsv donor
    Log.tsv tsv
    tsv.attach(donor, :fields => PCAWG::DONOR_FIELD).reorder(PCAWG::DONOR_FIELD, ["tumor_wgs_aliquot_id"]).to_single
  end

  PCAWG.claim PCAWG.donor_samples, :proc do
    tsv = PCAWG.donor_sample_info.tsv :fields => [SPECIMEN_FIELD, SAMPLE_FIELD, SNV_SAMPLE_FIELD_ORIG, EXPRESSION_SAMPLE_FIELD, NORMAL_EXPRESSION_SAMPLE_FIELD]
    preferred_samples = PCAWG.preferred_samples.tsv
    tsv.add_field SNV_SAMPLE_FIELD do |donor,values|
      ids = values[SNV_SAMPLE_FIELD_ORIG]
       if ids.length == 1
         ids
       else
         pref = preferred_samples[donor]
         good = (ids & [pref]).first
         iii [donor,ids,pref,good] if good.nil?
         [good]
       end
    end
    tsv.to_s
  end

  PCAWG.claim PCAWG.specimen_histology, :proc do
    histology = PCAWG::DATA_DIR["pcawg_specimen_histology_August2016_v1.tsv"].tsv :type => :double, :key_field => 'icgc_specimen_id', :sep2 => ',', :header_hash => '# '
    histology.to_s
  end

  PCAWG.claim PCAWG.subtype_info, :proc do
    tsv = PCAWG::DATA_DIR["tumour_subtype_consolidation_map.tsv - Unique List of Tumour Types_August.tsv"].tsv :type => :list, :key_field => "Abbreviation", :header_hash => ''
    tsv.to_s
  end

  def self.all_abbreviations
    PCAWG.specimen_histology.tsv(:key_field => "histology_abbreviation", :fields => [], :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache).keys
  end

  def self.all_histologies(tier = 1)
    PCAWG.specimen_histology.tsv(:key_field => "histology_tier"+tier.to_s, :fields => [], :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache).keys
  end

  def self.specimens_with_histology(h, key = nil)
    return all_specimens if h == "PCAWG"
    key, h = h.split ":" if h.include?(":") and key.nil?
    key ||= 'histology_abbreviation'
    PCAWG.specimen_histology.tsv.select(key => h, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache).keys
  end

  def self.all_specimens
    Persist.persist('all_specimens', :array, :dir => PCAWG::PROJECT_VAR_DIR.cache, :persist => true) do
      PCAWG.specimen_histology.tsv(:fields => [], :type => :single, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache).keys
    end
  end

  def self.specimen_histology_key(s, key = 'histology_abbreviation')
    @@speciment_histology_index ||= {}
    @@speciment_histology_index[key] ||= PCAWG.specimen_histology.tsv(:fields => [key], :key_field => "icgc_specimen_id", :type => :single, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache)
    @@speciment_histology_index[key][s]
  end

  def self.specimen_donor(specimen)
    @@index ||= PCAWG.donor_samples.index :target => DONOR_FIELD, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache, :data_persist => true
    donor = @@index[specimen]
    raise "No donor for specimen: #{ specimen }" if donor.nil?
    donor
  end

  def self.sample_donor(sample)
    @@index ||= PCAWG.donor_samples.index :target => DONOR_FIELD, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache
    @@index[sample]
  end

  def self.donors_with_histology(*args)
    specimens = specimens_with_histology(*args)
    donors = specimens.collect{|specimen| PCAWG.specimen_donor(specimen)}.uniq
    @@blacklisted_donors ||= PCAWG.blacklisted_donors.read.split("\n")
    donors = donors - @@blacklisted_donors
    Donor.setup(donors)
    donors.extend AnnotatedArray
    donors
  end

  def self.donor_histology(d, *args)
    specimens = Donor.donor2tumor_specimen[d]
    return specimens.collect{|s| specimen_histology_key(s, *args)}.compact.first
  end

  def self.abb_color(abb)
    @@color_index ||= PCAWG.subtype_info.tsv :key_field => 'Abbreviation', :fields => ["Color (RGB code)"], :type => :single
    @@color_index[abb]
  end
end
