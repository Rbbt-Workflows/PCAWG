Workflow.require_workflow "Study"
module Study

  property :donors => :single do
    samples = PCAWG.donors_with_histology(self)
    Sample.setup(samples, :cohort => self)
  end
  
  helper :organism do
    PCAWG.organism
  end

  task :organism => :string do
    PCAWG.organism
  end

  task :mappable_genes => :array do
    []
  end

  dep Sample, :genes_with_enhancer_mutations, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect{|s| s.genes_with_enhancer_mutations(:job, options) }
  end
  task :enhancer_mutation_incidence => :tsv do
    counts = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Sample"], :type => :flat)
    dependencies.each do |dep|
      sample = dep.clean_name
      TSV.traverse dep, :type => :array do |gene|
        counts[gene] ||= []
        counts[gene] << sample
      end
    end
    counts
  end

  dep Sample, :SV_candidate_genes, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect{|s| s.SV_candidate_genes(:job) }
  end
  task :SV_candidate_gene_incidence => :tsv do
    counts = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Sample"], :type => :flat)
    dependencies.each do |dep|
      sample = dep.clean_name
      TSV.traverse dep, :type => :array do |gene|
        counts[gene] ||= []
        counts[gene] << sample
      end
    end
    counts
  end

  dep Sample, :SV_candidate_fusions, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect{|s| s.SV_candidate_fusions(:job, options) }
  end

  task :SV_candidate_fusion_incidence => :tsv do
    counts = TSV.setup({}, :key_field => "Ensembl Gene ID - Ensembl Gene ID", :fields => ["Sample"], :type => :flat)
    dependencies.each do |dep|
      sample = dep.clean_name.split(":").last
      TSV.traverse dep, :type => :array do |gene|
        counts[gene] ||= []
        counts[gene] << sample
      end
    end
    counts
  end

  task :SV_summary => :tsv do
    require 'rbbt/matrix/barcode'
    require 'rbbt/statistics/fisher'

    matrix = study.matrix :gene_expression

    fields = []
    fields << "Gene 1 (Associated Gene Name)"
    fields << "Gene 2 (Associated Gene Name)"
    fields << "Fusion donors"
    fields << "Expression for Gene 1 in fusion donors"
    fields << "Expression for Gene 1 in rest of donors"
    fields << "Expression for Gene 2 in fusion donors"
    fields << "Expression for Gene 2 in rest of donors"
    fields << "Barcode for Gene 1 level in fusion donors"
    fields << "Barcode for Gene 1 level in rest of donors"
    fields << "Barcode for Gene 2 level in fusion donors"
    fields << "Barcode for Gene 2 level in rest of donors"
    fields << "Fisher p-value for gene 1"
    fields << "Fisher p-value for gene 2"
    fields << "Best Fisher p-value"


    tsv = TSV.setup({}, :key_field => "Candidate fusion", :fields => fields, :namespace => organism, :type => :double)

    if not matrix.nil? and matrix.samples and matrix.samples.any?
      all_rna_samples = matrix.samples

      barcode =  matrix.to_barcode_ruby.tsv 
      mat_values = matrix.tsv

      ensg2name = Organism.identifiers(organism).index :target => "Associated Gene Name", :persist => true
      donor2sample = PCAWG.donor_rna_samples.index :target => PCAWG::RNA_TUMOR_SAMPLE

      TSV.traverse study.SV_candidate_fusion_incidence, :bar => true, :into => tsv, :cpus => 10 do |fusion,donors|
        fusion = fusion.first if Array === fusion
        ens1, ens2 = fusion.split"-"
        next if ens1 == ens2
        next if donors.length < 3

        donors = donors.collect{|d| d.split(":").last }
        Sample.setup(donors, :cohort => study)

        samples = donor2sample.values_at *donors
        rna_seq_samples = all_rna_samples & samples

        barcode1 = barcode[ens1]
        barcode2 = barcode[ens2]
        next if barcode1.nil? or barcode2.nil?

        ens1_barcode_true = barcode1.values_at(*rna_seq_samples)
        ens1_barcode_false = barcode1.values_at(*(all_rna_samples - rna_seq_samples))
        ens2_barcode_true = barcode2.values_at(*rna_seq_samples)
        ens2_barcode_false = barcode2.values_at(*(all_rna_samples - rna_seq_samples))

        fusion_classes = all_rna_samples.collect{|s| samples.include?(s) ? 1 : 0 }
        expression_classes1 = all_rna_samples.collect{|s| barcode1[s] }
        expression_classes2 = all_rna_samples.collect{|s| barcode2[s] }
        pvalue1 = Fisher.test_classification(fusion_classes, expression_classes1)
        pvalue2 = Fisher.test_classification(fusion_classes, expression_classes2)

        name1, name2 = ensg2name.values_at ens1, ens2

        values = []
        values << [name1]
        values << [name2]
        values << donors
        values << mat_values[ens1].values_at(*(rna_seq_samples))
        values << mat_values[ens1].values_at(*(all_rna_samples - rna_seq_samples))
        values << mat_values[ens2].values_at(*(rna_seq_samples))
        values << mat_values[ens2].values_at(*(all_rna_samples - rna_seq_samples))
        values << barcode[ens1].values_at(*(rna_seq_samples))
        values << barcode[ens1].values_at(*(all_rna_samples - rna_seq_samples))
        values << barcode[ens2].values_at(*(rna_seq_samples))
        values << barcode[ens2].values_at(*(all_rna_samples - rna_seq_samples))
        values << [pvalue1]
        values << [pvalue2]
        values << [[pvalue2, pvalue1].min]

        [fusion, values]
      end
    end
    tsv
  end
  task :expression_samples => :array do
    donors = study.donors
    samples = Sample.setup(donors.expression_samples.flatten.compact, :cohort => self)
    samples.extend AnnotatedArray
    samples
  end


  property :organism => :single do
    PCAWG.organism
  end

  property :samples => :single do
    donors
  end

  property :genotyped_samples => :single do
    samples
  end

  property :cnv_samples => :single do
    donors.select(:has_cnv?)
  end

  property :sv_samples => :single do
    donors.select(:has_sv?)
  end


  property :has_cnv? => :single do
    cnv_samples.any?
  end

  property :has_sv? => :single do
    sv_samples.any?
  end

  property :has_drivers? => :single do
    !! PCAWG.driver_dir(self)
  end

  property :drivers => :single do |type,threshold=nil,fdr=nil|
    threshold = 0.1 if threshold.nil?
    fdr = true if fdr.nil?
    if has_drivers?
      tsv = PCAWG.driver_dir(self)[type].tsv :type => :list, :cast => :to_f
      tsv = FDR.adjust_hash!(tsv) if fdr
      if type == "enhancer"
        drivers = tsv.select("p-value"){|p| !p.nil? && (p.to_f < threshold.to_f) }.keys.collect{|g| g.split("::")[1].sub('-',':') }
        ChromosomeRange.setup(drivers, PCAWG.organism)
      else
        drivers = tsv.select("p-value"){|p| !p.nil? && (p.to_f < threshold.to_f) }.keys.collect{|g| g.split("::").last.gsub(/\.\d+$/,'') }
        Gene.setup(drivers, "Ensembl Gene ID", PCAWG.organism)
      end
    else
      nil
    end
  end

  property :has_candidate_drivers? => :single do
    !! PCAWG.candidate_driver_dir(self)
  end

  property :candidate_drivers => :single do |type|
    if has_candidate_drivers? and PCAWG.candidate_driver_dir(self)[type].exists?
      tsv = PCAWG.candidate_driver_dir(self)[type].tsv :type => :list, :fields => []
      tsv.keys.ensembl
    else
      []
    end
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
                                                 tsv.annotate data
                                                 data.filename = "Sample info #{study}"
                                                 data
                                               end
                                               info.entity_options = {:cohort => study}
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
    return nil if file.nil?
    sample_info = sample_info(study)
    value_type = study_info(study)[:expression_type]
    organism = study_info(study)[:organism]
    RbbtMatrix.new file.find, nil, value_type, format, organism
  end

  property :sample_extended_info => :single do
    donors = self.donors
    tsv = PCAWG.donor_clinical.tsv 
    tsv.select(PCAWG::DONOR_FIELD => donors)
  end

  dep Sample, :gene_extra_status, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    if study.has_cnv?
      study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study) unless Sample === sample; sample.gene_extra_status(:job, options) }.flatten
    else
      []
    end
  end
  task :sample_gene_extra => :tsv do
    if dependencies.any?
      parser = TSV::Parser.new dependencies.first
      fields = parser.fields
      fields.unshift "Sample"
      header = TSV.header_lines(parser.key_field, parser.fields, parser.options.merge(:type => :double))

      io = Misc.open_pipe do |sin|
        sin.puts header

        TSV.traverse dependencies, :type => :array do |job|
          sample = job.clean_name.split(":").last
          TSV.traverse job, :type => :array do |line|
            next if line =~ /^#/
            gene,_sep,status = line.partition("\t")
            parts = [gene, sample, status]
            sin.puts parts * "\t"
          end
        end
      end

      TSV.collapse_stream io
    else
      ""
    end
  end

  dep Sample, :gene_cnv_status_focal, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    if study.has_cnv?
      study.cnv_samples.collect{|sample| Sample.setup(sample, :cohort => study) unless Sample === sample; sample.gene_cnv_status_focal(:job, options) }.flatten
    else
      []
    end
  end
  task :sample_gene_cnvs_focal => :tsv do
    if dependencies.any?
      parser = TSV::Parser.new dependencies.first
      fields = parser.fields
      fields.unshift "Sample"
      header = TSV.header_lines(parser.key_field, parser.fields, parser.options.merge(:type => :double))

      io = Misc.open_pipe do |sin|
        sin.puts header

        TSV.traverse dependencies, :type => :array do |job|
          sample = job.clean_name.split(":").last
          TSV.traverse job, :type => :array do |line|
            next if line =~ /^#/
            gene,_sep,status = line.partition("\t")
            parts = [gene, sample, status]
            sin.puts parts * "\t"
          end
        end
      end

      TSV.collapse_stream io
    else
      ""
    end
  end

  dep :sample_gene_mutations, :compute => :produce
  dep :sample_gene_cnvs_focal, :compute => :produce
  task :sample_genes => :tsv do
    if study.has_cnv?
      io = TSV.paste_streams [step(:sample_gene_mutations), step(:sample_gene_cnvs_focal)]
      parser = TSV::Parser.new io
      dumper = TSV::Dumper.new parser.options.merge(:fields => parser.fields[0..-3] + parser.fields[-1..-1])
      dumper.init
      TSV.traverse parser, :into => dumper do |gene,values|
        gene = gene.first if Array === gene

        samples, *rest = values
        cnv = rest.pop
        cnv_samples = rest.pop
        new_values = rest
        new_cnv = ['normal'] * samples.length
        cnv_samples.each_with_index do |cnv_sample,i|
          index = samples.index cnv_sample
          if index.nil?
            samples << cnv_sample
            new_values.each{|l| l << 'false'}
            new_cnv << cnv[i]
          else
            new_cnv[index] = cnv[i]
          end
        end
        new_values.unshift(samples)
        new_values.push(new_cnv)
        [gene, new_values]
      end
    else
      TSV.get_stream step(:sample_gene_mutations)
    end
  end
  
end

Study.update_task_properties
