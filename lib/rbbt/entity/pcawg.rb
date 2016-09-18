require 'rbbt/entity/pcawg/study'
require 'rbbt/entity/pcawg/sample'

module PCAWG

  def self.all_donors
    @@all_donors ||= PCAWG.donor_wgs_samples.tsv(:key_field => PCAWG::DONOR_FIELD, :fields => [], :type => :list).keys
  end

  def self.donors_with_histology(s)
    if s.include? "="
      key, value = s.split("=")
    else
      key, value = 'histology_abbreviation', s
    end

    @@good_donors ||= PCAWG.donor_wgs_samples.tsv.keys

    if s == "PCAWG"
      @@good_donors
    else
      @@donor_histology ||= PCAWG.donor_histology.tsv 
      donors = @@good_donors
      donors &= @@donor_histology.select(key => value).keys 
      donors
    end
  end

  def self.all_abbreviations
    PCAWG.donor_histology.tsv(:key_field => "histology_abbreviation", :fields => [], :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache).keys
  end

  def self.all_histologies(tier = 1)
    PCAWG.donor_histology.tsv(:key_field => "histology_tier"+tier.to_s, :fields => [], :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache).keys
  end

end
