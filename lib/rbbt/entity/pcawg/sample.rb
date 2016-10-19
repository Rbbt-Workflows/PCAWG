Workflow.require_workflow "Sample"

module Sample

  helper :watson do
    true
  end

  helper :organism do
    PCAWG.organism
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

    dumper = TSV::Dumper.new(:key_field => "SV", :fields => ["Ensembl Gene ID before left boundary", "Ensembl Gene ID after left boundary","Ensembl Gene ID before right boundary","Ensembl Gene ID after right boundary"], :organism => organism, :type => :double)
    dumper.init
    TSV.traverse Sequence.job(:genes_at_ranges, clean_name, :ranges => stream, :organism => organism).produce(true), :into => dumper do |range, genes|
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
      line.split("\t").collect{|v| v.split("|").reject{|v| v.empty?} * "|" } * "\t"
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
      next if left.empty? or right.empty?
      pairs = left.collect{|l| right.collect{|r| l == r ? nil : [l,r].sort * "-" } }.flatten.compact
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

  property :expression_samples => :array do
    index = PCAWG.donor_rna_samples.index :target => PCAWG::RNA_TUMOR_SAMPLE
    samples = index.chunked_values_at(self).compact
    Sample.setup(samples, :cohort => cohort)
  end

  property :abbr_color => :single do
    raise "TODO"
    abbr = PCAWG.donor_histology(self, 'histology_abbreviation')
    PCAWG.abb_color(abbr)
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
    false
  end
end

Sample.update_task_properties
