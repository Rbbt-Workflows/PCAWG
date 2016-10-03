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
