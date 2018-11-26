module PCAWG

  WGS_SAMPLE_FIELD = 'tumor_wgs_aliquot_id'
  DONOR_FIELD = 'icgc_donor_id'
  RNA_TUMOR_SAMPLE = 'RNA_tumor_sample'
  RNA_NORMAL_SAMPLE = 'RNA_normal_sample'

  PCAWG.claim PCAWG.donor_specimens, :proc do |filename|
    tsv = DATA_DIR['release_may2016.v1.4.tsv'].tsv :key_field => DONOR_FIELD, :fields => %w(normal_wgs_icgc_specimen_id tumor_wgs_icgc_specimen_id normal_rna_seq_icgc_specimen_id tumor_rna_seq_icgc_specimen_id), :header_hash => '', :sep2 => ','
    tsv.fields = ["Normal WGS Specimen", "Tumor WGS Specimen", "Normal RNA Specimen", "Tumor RNA Specimen"]
    tsv.to_s
  end

  PCAWG.claim PCAWG.selected_donor_samples, :proc do |filename|
    list = TSV.traverse DATA_DIR["aliquot_donor_tumor.whitelist.tsv.gz"].find, :into => [], :type => :array do |line|
      next if line =~ /DONOR_UNIQ/
      line.split("\t")[1]
    end
    list * "\n" + "\n"
  end

  PCAWG.claim PCAWG.donor_wgs_samples, :proc do |filename|
    tsv = DATA_DIR['release_may2016.v1.4.tsv'].tsv :key_field => WGS_SAMPLE_FIELD, :fields => [DONOR_FIELD], :header_hash => '', :sep2 => ','
    selected_donor_samples = PCAWG.selected_donor_samples.list
    tsv = tsv.select selected_donor_samples
    tsv = tsv.reorder 'icgc_donor_id', %w(tumor_wgs_aliquot_id)
    tsv.to_single.to_s
  end
  
  PCAWG.claim PCAWG.donor_rna_samples, :proc do |filename|
    tsv = DATA_DIR['release_may2016.v1.4.tsv'].tsv :key_field => DONOR_FIELD, :fields => %w(tumor_rna_seq_aliquot_id normal_rna_seq_aliquot_id), :header_hash => '', :sep2 => ','

    tsv.fields = [RNA_TUMOR_SAMPLE, RNA_NORMAL_SAMPLE]

    tsv.to_list.to_s
  end

  PCAWG.claim PCAWG.donor_clinical, :proc do |filename|
    fields = %w(donor_sex donor_vital_status donor_diagnosis_icd10 first_therapy_type first_therapy_response donor_age_at_diagnosis donor_survival_time donor_interval_of_last_followup tobacco_smoking_history_indicator tobacco_smoking_intensity alcohol_history alcohol_history_intensity)
    tsv = DATA_DIR["pcawg_donor_clinical_August2016_v9.tsv"].tsv :key_field => 'icgc_donor_id', :fields => fields, :type => :list, :fix => Proc.new{|l| l.gsub(/([a-z]+)CA(\s)/, '\1Ca\2')}
    tsv.to_s
  end

  PCAWG.claim PCAWG.donor_ICGC_therapy, :proc do 
    name = "donor_therapy.all_projects.tsv.gz"
    file = "/current/Summary/" << name
    url = "http://dcc.icgc.org/api/v1/download?fn=" << file
    tsv = TSV.open(url, :header_hash => "")
    good_fields = tsv.fields - %w(project_code submitted_donor_id)
    all_donors = PCAWG.donor_specimens.tsv.keys
    tsv.slice(good_fields).select(all_donors).to_s
  end

  PCAWG.claim PCAWG.donor_ICGC_surgery, :proc do 
    name = "donor_surgery.all_projects.tsv.gz"
    file = "/current/Summary/" << name
    url = "http://dcc.icgc.org/api/v1/download?fn=" << file
    tsv = TSV.open(url, :header_hash => "")
    good_fields = tsv.fields - %w(project_code submitted_donor_id)
    all_donors = PCAWG.donor_specimens.tsv.keys
    tsv.slice(good_fields).select(all_donors).to_s
  end

  PCAWG.claim PCAWG.donor_ICGC_exposure, :proc do 
    name = "donor_exposure.all_projects.tsv.gz"
    file = "/current/Summary/" << name
    url = "http://dcc.icgc.org/api/v1/download?fn=" << file
    tsv = TSV.open(url, :header_hash => "")
    good_fields = tsv.fields - %w(project_code submitted_donor_id)
    all_donors = PCAWG.donor_specimens.tsv.keys
    tsv.slice(good_fields).select(all_donors).to_s
  end

  PCAWG.claim PCAWG.donor_ICGC_family, :proc do 
    name = "donor_family.all_projects.tsv.gz"
    file = "/current/Summary/" << name
    url = "http://dcc.icgc.org/api/v1/download?fn=" << file
    tsv = TSV.open(url, :header_hash => "")
    good_fields = tsv.fields - %w(project_code submitted_donor_id)
    all_donors = PCAWG.donor_specimens.tsv.keys
    tsv.slice(good_fields).select(all_donors).to_s
  end

  PCAWG.claim PCAWG.donor_histology, :proc do |filename|
    fields = %w(organ_system histology_abbreviation histology_tier1 histology_tier2 histology_tier3 histology_tier4 tumour_histological_code tumour_histological_type tumour_stage tumour_grade specimen_donor_treatment_type)
    tsv = DATA_DIR['pcawg_specimen_histology_August2016_v9.tsv'].tsv :key_field => 'icgc_donor_id', :fields => fields, :type => :double, :sep2 => /,\s*/, :fix => Proc.new {|l| l.gsub(/([a-z]+)CA(\s)/, '\1Ca\2')}
  end

  PCAWG.claim PCAWG.donor_meta_cohort, :proc do |filename|
    sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
    donor_groups = TSV.setup({}, :key_field => PCAWG::DONOR_FIELD, :fields => ["Meta Cohort"], :type => :flat)
    DATA_DIR.meta_tumor_cohorts.glob("*.whitelist.*").each do |file|
      name = File.basename(file).sub(".whitelist.txt",'')
      samples = Open.read(file).split("\n").collect{|l| l.split("\t").first}
      donors = sample2donor.values_at *samples
      donors.compact.each do |donor|
        donor_groups[donor] ||= []
        donor_groups[donor] << name
      end
    end
    donor_groups.to_s
  end

  PCAWG.claim PCAWG.donor_other_cohorts["PANCANCER_no_melanoma_lymph"], :proc do 
    all_donors = PCAWG.all_donors.sort
    lymph = PCAWG.donors_with_histology("meta=Lymph_tumors")
    mela = PCAWG.donors_with_histology("Skin-Melanoma")
    all_donors - lymph - mela
  end

  PCAWG.claim PCAWG.donor_other_cohorts["PANCANCER_no_melanoma_lymph.NOHYPER"], :proc do 
    PCAWG.donor_other_cohorts["PANCANCER_no_melanoma_lymph"].list & PCAWG.all_donors(true)
  end

  def self.all_donors(remove_hyper = false, hyper_threshold = 90_000)
    @@all_donors ||= PCAWG.donor_wgs_samples.tsv(:key_field => PCAWG::DONOR_FIELD, :fields => [], :type => :list).keys
    if remove_hyper
      Sample.setup(@@all_donors, :cohort => "PCAWG").select{|s| s.num_genomic_mutations < hyper_threshold}
    else
      @@all_donors
    end
  end

  def self.donors_with_histology(s)
    orig = s
    if s =~ /(.*).NOHYPER$/
      remove_hyper = true
      s = $1
    end

    if s.include? "="
      key, value = s.split("=")
    else
      key, value = 'histology_abbreviation', s
    end

    all_donors = PCAWG.all_donors(remove_hyper)

    if s == "PCAWG"
      all_donors
    elsif key == "meta"
      PCAWG.donor_meta_cohort.tsv.select("Meta Cohort" => value).keys & all_donors
    elsif key == "other"
      PCAWG.donor_other_cohorts[value].list & all_donors
    else
      begin
        donors = PCAWG.donor_histology.tsv.select(key => value).keys
        raise "No donors" if donors.empty?
        donors & all_donors
      rescue
        PCAWG.donor_other_cohorts[s].list & all_donors
      end
    end
  end

  def self.all_abbreviations
    PCAWG.donor_histology.tsv(:key_field => "histology_abbreviation", :fields => [], :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache).keys
  end

  def self.all_meta
    PCAWG.donor_meta_cohort.tsv.values.flatten.compact.uniq
  end

  def self.all_histologies(tier = 1)
    PCAWG.donor_histology.tsv(:key_field => "histology_tier"+tier.to_s, :fields => [], :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache).keys
  end

end
