
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
    Donor.setup(PCAWG.donors_with_histology(self), :cohort => "PCAWG")
  end

  property :genotyped_samples => :single do
    samples = Sample.setup(donors.collect{|donor| donor.SNV_sample}.compact.flatten)
    samples.extend AnnotatedArray
    samples.cohort = "PCAWG"
    samples
  end

  property :expression_samples => :single do
    samples = Sample.setup(donors.collect{|donor| donor.expression_samples}.compact.flatten)
    samples.extend AnnotatedArray
    samples.cohort = "PCAWG"
    samples
  end

  property :has_cnv? => :single do
    false
  end

  property :condition => :single do
    self
  end

  property :metadata => :single do
    {}
  end

  def self.study_info(study)
    @@study_info ||= {:watson => true, :organism => PCAWG.organism, :expression_type => "RPKM"}
  end

  def self.main_donor_info
    Persist.persist_tsv nil, 'PCAWG', {}, :dir => PCAWG::PROJECT_VAR_DIR.donor_info, :persist => true do |data|
      tsv = PCAWG.donor_samples.tsv
      tsv = tsv.attach Rbbt.data.preliminary_final_release["pcawg_donor_clinical_May2016_v2.tsv"].find
      tsv = tsv.attach Rbbt.data.preliminary_final_release["pcawg_specimen_histology_May2016_v2.tsv"].find
      ppp tsv.fields * "\n"
      data.serializer = :double
      data.merge! tsv
      tsv.annotate data
      data
    end
  end

  def self.donor_info(study)
    return main_donor_info if study == "PCAWG"
    @@donor_info ||= {}
    @@donor_info[study] ||= begin
                              info = Persist.persist_tsv nil, study, {}, :dir => PCAWG::PROJECT_VAR_DIR.donor_info, :persist => true do |data|
                                study = Study.setup(study.dup)
                                tsv = main_donor_info.select study.donors
                                data.serializer = :double
                                tsv.annotate data
                                data.merge! tsv
                                data.filename = "Donor info #{study}"
                                data
                              end
                              info.entity_options = {:cohort => "PCAWG"}
                              info
                            end
  end

  def self.sample_info(study, sample_field = PCAWG::SNV_SAMPLE_FIELD)
    @@sample_info ||= {}
    @@sample_info[[study, sample_field]] ||= begin 
                                               info = Persist.persist_tsv nil, study, {}, :dir => PCAWG::PROJECT_VAR_DIR.sample_info, :persist => true do |data|
                                                 tsv = donor_info(study)
                                                 tsv = tsv.select(study.donors) unless study == "PCAWG"
                                                 tsv = tsv.reorder sample_field
                                                 data.serializer = :double
                                                 data.merge! tsv
                                                 data.filename = "Sample info #{study}"
                                                 tsv.annotate data
                                                 data
                                               end
                                               info.entity_options = {:cohort => "PCAWG"}
                                               info
                                             end
  end

  def self.matrix_file(study, matrix)
    if study.expression_samples.any?
      file = PCAWG::PROJECT_VAR_DIR.matrices[study].find
      Persist.persist study, :tsv, :file => file, :dir => PCAWG::PROJECT_VAR_DIR.matrices, :persist => true, :no_load => true do 
        orig_file = Rbbt.data.preliminary_final_release["joint_fpkm_uq.tsv"].find
        fields = TSV.parse_header(orig_file, :header_hash => '').fields
        study = Study.setup(study.dup)
        pos = study.expression_samples.collect{|s| fields.index s}.compact.collect{|i| i + 2}

        cmd = "zcat #{orig_file} | cut -f 1,#{pos * ","} |sed 's/\\(ENSG[[:digit:]]\\+\\)\\.[[:digit:]]\\+/\\1/' "
        tsv = TSV.open(CMD.cmd(cmd, :pipe => true), :header_hash => '')
        tsv.key_field = "Ensembl Gene ID"
        tsv.to_s
      end
      file
    else
      nil
    end
  end

  def self.matrices(study)
    if matrix_file(study, :gene_expression)
      [:gene_expression]
    else
      []
    end
  end


  def self.matrix(study, matrix, format = nil)
    file = matrix_file(study, matrix)
    sample_info = sample_info(study)
    value_type = study_info(study)[:expression_type]
    organism = study_info(study)[:organism]
    Matrix.new file.find, nil, value_type, format, organism
  end

  property :sample_extended_info => :single do
    donors = self.donors
    tsv = PCAWG.donor_clinical.tsv 
    tsv.select(PCAWG::DONOR_FIELD => donors).attach(PCAWG.donor_samples, :fields => [PCAWG::SNV_SAMPLE_FIELD]).reorder(PCAWG::SNV_SAMPLE_FIELD)
  end
  
  property :cnv_samples => :single do
    donors.select(:has_cnv?)
  end

  task :test => :boolean do
    study = clean_name.dup
    study = Study.setup(study)
    iii study.expression_samples
    raise "STOP"
    Log.tsv Study.sample_info(study)
    study.genotyped_samples.each{|s|
      iif [s, s.donor.expression_samples]
    }
    raise "Stop"
  end
end

Study.update_task_properties
