Workflow.require_workflow "Sample"
Workflow.require_workflow "Sequence"

module Sample

  helper :watson do
    true
  end

  helper :organism do
    PCAWG.organism
  end

  helper :clean_donor do
    sample.split(":").last
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

  dep :genomic_mutations, :compute => :produce
  dep Sequence, :intersect_bed, :positions => :genomic_mutations
  task :intersect_bed => :tsv do
    TSV.get_stream step(:intersect_bed)
  end

  dep :genomic_mutations, :compute => :produce
  dep Sequence, :intersect_bed, :sorted => true, :positions => :genomic_mutations, :bed_file => 'placeholder' do |jobname, options|
    PCAWG.regions.produce.glob("*").collect do |file|
      inputs = options.merge({:bed_file => "" + file.find.to_s})
      {:workflow => Sequence, :task => :intersect_bed, :jobname => jobname, :inputs => inputs}
    end
  end
  task :bed_region_mutations => :tsv do
    ios = []
    names = []
    dependencies.each do |dep|
      next unless dep.task_name == :intersect_bed
      name = File.basename dep.recursive_inputs.to_hash[:bed_file]
      io = Misc.collapse_stream(TSV.get_stream(dep))
      ios << io
      names << name
    end
    
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => names, :type => :double, :namespace => PCAWG.namespace
    dumper.init 

    TSV.traverse TSV.paste_streams(ios, :sort => false), :type => :array, :into => dumper, :bar => "Merging intersections with BED files" do |line|
      next if line =~ /^#/
      key, *values = line.split("\t")
      
      [key, values]
    end
  end

  dep :bed_region_mutations
  task :gene_bed_intersections => :tsv do
    parser = TSV::Parser.new step(:bed_region_mutations)
    fields = parser.fields
    dumper = TSV::Dumper.new(:key_field => "Ensembl Gene ID", :fields => ["Region type"], :type => :flat, :namespace => PCAWG.organism)
    dumper.init
    TSV.traverse parser, :into => dumper do |mutation, values|
      res = []
      fields.zip(values).each do |field, list|
        next if list.nil?
        genes = list.collect{|e| e.split("::").select{|_e| _e =~ /ENSG/}.first }.compact
        genes = genes.collect{|e| e.split(".").first }
        genes.each do |gene|
          res << [gene, [field]]
        end
      end
      next if res.empty?
      res.extend MultipleResult
      res
    end

    io = TSV.collapse_stream dumper.stream
    TSV.traverse io, :into => :stream, :type => :array do |line|
      line.gsub('|', "\t")
    end
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

  property :get_genes_extra => :single do |type|
    genes = self.gene_extra_status.select(type => "true").keys
    Gene.setup(genes.dup, "Ensembl Gene ID", organism).extend AnnotatedArray
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

  dep :organism
  dep :gene_cnv_status
  task :gene_cnv_status_focal => :tsv do 
    index = PCAWG.donor_wgs_samples.index :target => PCAWG::WGS_SAMPLE_FIELD, :persist => true
    name2ensembl = Organism.identifiers(PCAWG.organism).index :target => "Ensembl Gene ID", :fields => "Associated Gene Name", :persist => true
    sample = index[clean_donor]
    begin
      focal = PCAWG::DATA_DIR["focal_data_by_genes.rmcnv.pt_170207.txt"].tsv(:type => :list, :cast => :to_f, :persist => true, :header_hash => '').select(sample){|v| v != 0 }.keys
    rescue
      focal = []
    end
    focal_ensembl = name2ensembl.chunked_values_at focal
    step(:gene_cnv_status).load.select(focal_ensembl)
  end

  input :gain_cnv_threshold, :float, "Copy number threshold to consider the gene", 1
  task :gained_GISTIC => :array do |threshold|
    name2ensg = Organism.identifiers(PCAWG.organism).index :target => "Ensembl Gene ID", :persist => true, :order => true
    donor = clean_donor
    genes = begin
      PCAWG.matrices.copy_number.tsv(:type => :list, :cast => :to_f, :persist => true).select(donor){|v| v > threshold}.keys
    rescue
      []
    end
    name2ensg.values_at(*genes).compact
  end

  input :loss_cnv_threshold, :float, "Copy number threshold to consider the gene", -1
  task :lost_GISTIC => :array do |threshold|
    name2ensg = Organism.identifiers(PCAWG.organism).index :target => "Ensembl Gene ID", :persist => true, :order => true
    donor = clean_donor
    genes = begin
      PCAWG.matrices.copy_number.tsv(:type => :list, :cast => :to_f, :persist => true).select(donor){|v| v < threshold}.keys
    rescue
      Log.exception $!
      []
    end
    name2ensg.values_at(*genes).compact
  end

  input :gain_cnv_threshold, :float, "Copy number threshold to consider the gene", 1
  task :gained_focal_GISTIC => :array do |threshold|
    name2ensg = Organism.identifiers(PCAWG.organism).index :target => "Ensembl Gene ID", :persist => true, :order => true
    donor = clean_donor
    genes = begin
      PCAWG.matrices.copy_number_focal.tsv(:type => :single, :cast => :to_f, :persist => true).select(donor){|v| v > threshold}.keys
    rescue
      []
    end
    name2ensg.values_at(*genes).compact
  end

  input :loss_cnv_threshold, :float, "Copy number threshold to consider the gene", -1
  task :lost_focal_GISTIC => :array do |threshold|
    name2ensg = Organism.identifiers(PCAWG.organism).index :target => "Ensembl Gene ID", :persist => true, :order => true
    donor = clean_donor
    genes = begin
      PCAWG.matrices.copy_number_focal.tsv(:type => :single, :cast => :to_f, :persist => true).select(donor){|v| v < threshold}.keys
    rescue
      []
    end
    name2ensg.values_at(*genes).compact
  end

  property :pcawg_mutation_signature do
    donor = self
    tsv = PCAWG.mutation_signature.assignments.tsv
    raise "No signature for donor #{ donor }" unless tsv.include?(donor)
    tsv.select(donor)
  end

  dep :lost_focal_GISTIC, :compute => :produce
  dep :gained_focal_GISTIC, :compute => :produce
  dep :bed_region_mutations
  task :gene_extra_status => :tsv do
    lost = step(:lost_focal_GISTIC).load
    gained = step(:gained_focal_GISTIC).load

    parser = TSV::Parser.new step(:bed_region_mutations)
    fields = parser.fields
    all_fields = fields + ["Focal GISTIC Lost", "Focal GISTIC Gained"]

    all_fields = all_fields.reject{|f| f.include? 'enhanc' }
    all_fields = all_fields.reject{|f| f == 'gc19_pc.cds.bed' }

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => all_fields, :type => :double, :namespace => PCAWG.organism)

    lost_pos = all_fields.index "Focal GISTIC Lost"
    lost.each do |gene|
      val = ['false'] * all_fields.length
      val[lost_pos] = 'true'
      tsv.zip_new gene, val
    end

    gained_pos = all_fields.index "Focal GISTIC Gained"
    gained.each do |gene|
      val = ['false'] * all_fields.length
      val[gained_pos] = 'true'
      tsv.zip_new gene, val
    end


    TSV.traverse parser, :into => tsv do |mutation,values_list|
      gene_matches = {}
      fields.zip(values_list).each do |field, values|
        next if values.nil? or values.empty?
        genes = values.collect{|v| v.split("::").select{|v| v =~ /ENSG/}}.flatten.uniq
        genes.each do |gene|
          gene.sub!(/\.\d+/,'')
          gene_matches[gene] ||= []
          gene_matches[gene] << field
        end
      end
      res = []
      gene_matches.each do |gene,matches|
        r = all_fields.collect do |field|
          matches.include?(field) ? 'true' : 'false'
        end
        res << [gene, r]
      end
      res.extend MultipleResult

      res
    end

    dumper2 = TSV::Dumper.new :key_field => "Ensembl Gene ID", :fields => all_fields, :type => :list, :namespace => PCAWG.organism
    dumper2.init
    TSV.traverse tsv, :type => :array, :into => dumper2 do |gene, values|
      values = values.collect{|l| 
        l.include?('true') ? 'true' : 'false' 
      } 
      [gene, values]
    end
  end

  dep :organism
  dep :watson
  dep :genomic_mutations, :compute => :produce
  dep Sequence, :affected_genes, :mutations => :genomic_mutations, :organism => :organism, :watson => :watson
  task :gene_timing => :tsv do
    organism = step(:organism).load
    timing = PCAWG.clonality.timing.produce[clean_donor].tsv :unnamed => true

    gene_timings = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Timing"], :type => :flat, :namespace => organism)
    TSV.traverse step(:affected_genes) do |mutation, genes|
      mutation = mutation.first if Array === mutation
      chr, pos = mutation.split(":").values_at(0,1)
      time = timing[[chr, pos] * ":"]
      time = timing[[chr, (pos.to_i + 1).to_s] * ":"] if time.nil?
      genes.each do |gene|
        gene_timings[gene] ||= []
        gene_timings[gene] << time
      end
    end
    gene_timings
  end
end

Sample.update_task_properties
