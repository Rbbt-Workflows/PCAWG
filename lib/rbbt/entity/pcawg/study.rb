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
    Sample.setup(PCAWG.donors_with_histology(self), :cohort => "PCAWG")
  end
  
  property :samples => :single do
    donors
  end

  property :genotyped_samples => :single do
    donors
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
    @@study_info ||= {:watson => true, :organism => PCAWG.organism, :expression_type => "FPKM"}
  end

  def self.sample_info(study, sample_field = PCAWG::DONOR_FIELD)
    @@sample_info ||= {}
    @@sample_info[[study, sample_field]] ||= begin 
                                               info = Persist.persist_tsv nil, study, {}, :dir => PCAWG::PROJECT_VAR_DIR.sample_info, :persist => true do |data|
                                                 Study.setup(study) unless Study === study
                                                 tsv = PCAWG.donor_wgs_samples.tsv :type => :double
                                                 tsv = tsv.select(study.donors) unless study == "PCAWG"
                                                 tsv = tsv.attach PCAWG.donor_rna_samples
                                                 tsv = tsv.attach PCAWG.donor_histology
                                                 tsv = tsv.attach PCAWG.donor_clinical
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
        orig_file = PCAWG.matrices.gene_expression.find
        fields = TSV.parse_header(orig_file, :header_hash => '').fields
        study = Study.setup(study.dup)
        pos = study.expression_samples.collect{|s| fields.index s}.compact.collect{|i| i + 2}

        cmd = "cat #{orig_file} | cut -f 1,#{pos * ","} |sed 's/\\(ENSG[[:digit:]]\\+\\)\\.[[:digit:]]\\+/\\1/' "
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
    tsv.select(PCAWG::DONOR_FIELD => donors)
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
