Workflow.require_workflow "Sample"

module Sample

  helper :watson do
    true
  end

  helper :organism do
    PCAWG.organism
  end

  helper :cnv_status do |info|
    total, mj, mn, stars = info
    statuses = []
    case total.to_i
    when 0
      'complete_loss'
    when 1
      'loss'
    when 2
      if mj == '2'
        'LOH'
      else
        nil
      end
    when 3
      'gain'
    else
      if total.to_i > 3
        'big_gain'
      elsif total.to_i < 0
        'complete_loss'
      else
        raise "Unknown CNV: #{ total }"
      end
    end
  end

  property :sample_code => :single do
    clean = self.split(":").last
    "PCAWG:" << clean
  end

  task :homozygous => :array do
    []
  end

  task :homozygous_genes => :array do
    []
  end

  task :extended_vcf => :tsv do
    nil
  end

  task :organism => :string do
    PCAWG.organism
  end

  task :genomic_mutations => :array do
    Misc.sort_mutation_stream PCAWG.genotypes.produce[clean_name.split(":").last].open
  end

  task :cnvs => :array do
    Misc.sort_mutation_stream PCAWG.CNV.produce[clean_name.split(":").last].open
  end

  dep :genomic_mutations
  dep Sequence, :intersect_bed, :positions => :genomic_mutations
  task :intersect_bed => :tsv do
    TSV.get_stream step(:intersect_bed)
  end

  dep :intersect_bed, :bed_file => PCAWG.enhancer_ranges.produce
  task :genes_with_enhancer_mutations => :array do
    TSV.traverse step(:intersect_bed), :type => :array, :into => :stream do |line|
      genes = line.split("\t")[1].split(":")[3].split(";")
      genes.extend MultipleResult
    end
  end

  returns "StructuralVariant"
  task :SV => :array do
    file = PCAWG.SV.produce[clean_name.split(":").last]
    if file.exists?
      PCAWG.SV.produce[clean_name.split(":").last].open
    else
      []
    end
  end

  dep :SV
  input :size, :integer, "Number of bases at either side of boundary", 1000
  task :SV_boundaries => :tsv do |size|
    stream = TSV.traverse step(:SV), :type => :array, :into => :stream do |line|
      chr1,pos1,type,chr2,pos2 = line.chomp.split(":")
      pos1 = pos1.to_i
      pos2 = pos2.to_i
      boundaries = []
      boundaries << [chr1,pos1 - size,pos1, "pre", pos1, line] * ":"
      boundaries << [chr1,pos1,pos1 + size, "post", pos1, line] * ":"
      boundaries << [chr2,pos2 - size,pos2, "pre", pos2, line] * ":"
      boundaries << [chr2,pos2,pos2 + size, "post", pos2, line] * ":"
      boundaries.extend MultipleResult
      boundaries
    end

    dumper = TSV::Dumper.new(:key_field => "Structural Variant", :fields => ["Gene (Ensembl Gene ID) before left boundary", "Gene (Ensembl Gene ID) after left boundary","Gene (Ensembl Gene ID) before right boundary","Gene (Ensembl Gene ID) after right boundary"], :namespace => organism, :type => :double)
    dumper.init
    TSV.traverse Sequence.job(:genes_at_ranges, clean_name, :ranges => stream, :namespace => organism).produce(true), :into => dumper do |range, genes|
      range = range.first if Array === range
      parts = range.split(":")
      chr1, left, right, loc, pos1, *line = parts
      dir = pos1 == line[1] ? "left" : "right"
      sv = line  * ":"
      case dir
      when "left"
        case loc
        when "pre"
          [sv, [genes, nil, nil, nil]]
        when "post"
          [sv, [nil, genes, nil, nil]]
        end
      when "right"
        case loc
        when "pre"
          [sv, [nil, nil, genes, nil]]
        when "post"
          [sv, [nil, nil, nil, genes]]
        end
      end
    end
  
    stream = TSV.collapse_stream(dumper, :zipped => false)
    TSV.traverse stream, :into => :stream, :type =>:array do |line|
      parts = line.split("\t").collect{|v| v.split("|").reject{|v| v.empty?} * "|" } 
      next if parts.length == 1
      parts * "\t"
    end
  end

  dep :SV_boundaries
  input :translate, :boolean, "Translate gene IDS into symbol", false
  task :SV_candidate_fusions => :array do |translate|
    index = Organism.identifiers(organism).index :target => "Associated Gene Name", :persist => true
    stream = TSV.traverse step(:SV_boundaries), :into => :stream do |sv,values|
      next if values.nil?
      left = (values[0..1] || []).flatten.compact.uniq
      left = index.values_at *left if translate
      right = (values[2..3] || []).flatten.compact.uniq
      right = index.values_at *right if translate
      common = left & right
      left -= common
      right -= common
      next if left.compact.empty? or right.compact.empty?
      pairs = left.collect{|l| right.collect{|r| l == r ? nil : [l,r].sort * "-" } }.flatten.compact
      next if pairs.empty?
      pairs * "\n"
    end
    CMD.cmd('sort -u', :in => stream, :pipe => true, :no_fail => true)
  end

  dep :SV_boundaries
  task :SV_candidate_genes => :array do 
    stream = TSV.traverse step(:SV_boundaries), :into => :stream do |sv,values|
      values.flatten.compact.uniq * "\n"
    end
    CMD.cmd('sort -u', :in => stream, :pipe => true, :no_fail => true)
  end

  input :limit_to_tumor_type, :boolean, "Consider only genes considered drivers in tumor type", true
  input :use_cnv, :boolean, "Consider gene in CNV as drivers", false
  task :driver_genes => :array do |limit,cnv|
    s = sample.split(":").last
    sample = Sample.setup(s)
    abbr = sample.abbr
    genes = limit ? PCAWG.abbr_driver_genes[abbr] : PCAWG.all_driver_genes

    raise RbbtException, "Abbreviation #{abbr} not in driver list - #{[abbr, sample, limit, cnv]}" if genes.nil?

    affected_genes = sample.get_genes(:affected)


    affected_drivers_abbr = affected_genes & genes

    if cnv 
      (affected_drivers_abbr + (sample.get_genes(:gained) & genes & PCAWG.activating_driver_genes) + (sample.get_genes(:lost) & genes & PCAWG.lof_driver_genes)).uniq
    else
      affected_drivers_abbr
    end
  end

  property :expression_samples => :array do
    index = PCAWG.donor_rna_samples.tsv :fields => [PCAWG::RNA_TUMOR_SAMPLE], :type => :flat
    self.collect do |donor|
      code = donor.split(":").last
      s = index[code]
      Sample.setup(s, :cohort => cohort)
      s
    end
  end

  property :abbr => :single do
    (@@donor_histology_abbr ||= PCAWG.donor_histology.tsv(:fields => ["histology_abbreviation"], :type => :single))[self]
  end

  property :abbr_color => :single do
    raise "TODO"
    PCAWG.abb_color(self.abbr)
  end

  property :gene_status => :single2array do |genes|
    raise "TODO"
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

  property :has_expression? => :single do
    ! self.expression_samples.nil?
  end

  property :has_genotype? => :single do
    true
  end

  property :has_cnv? => :single do
    PCAWG.CNV.produce[self].exists?
  end

  property :has_sv? => :single do
    PCAWG.SV.produce[self].exists?
  end
  
  task :pandrugs => :tsv do 
    assoc_fields = %w(affected_gene gscore alteration target_marker resistance relation gene_symbol source)
    associations = PCAWG.pandrugs[sample].tsv :type => :double, :merge => true, :sep2 => '--NONE--', :key_field => 'show_drug_name', :fields => assoc_fields

    ppp associations.to_s
    Log.tsv associations

    info_fields = %w(score status pathology cancer extra extra2)
    drug_info = PCAWG.pandrugs[sample].tsv :type => :double, :merge => true, :sep2 => '--NONE--', :key_field => 'show_drug_name', :fields => info_fields

    ppp drug_info.to_s
    Log.tsv drug_info

    raise 
  end

  input :gain_cnv_threshold, :float, "Copy number threshold to consider the gene", 1
  task :gained_GISTIC => :array do |threshold|
    donor = sample
    begin
      PCAWG.matrices.copy_number.tsv(:fields => [donor], :type => :single, :cast => :to_f).select(donor){|v| v > threshold}.keys
    rescue
      []
    end
  end

  input :loss_cnv_threshold, :float, "Copy number threshold to consider the gene", -1
  task :lost_GISTIC => :array do |threshold|
    donor = sample
    begin
      PCAWG.matrices.copy_number.tsv(:fields => [donor], :type => :single, :cast => :to_f).select(donor){|v| v < threshold}.keys
    rescue
      []
    end
  end
end

Sample.update_task_properties
