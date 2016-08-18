
module Donor
  extend Entity

  annotation :cohort

  property :histology_abbreviation => :single do
    specimens = PCAWG.donor2tumor_specimen[self]
    index = PCAWG.specimen_histology.tsv :key_field => 'icgc_specimen_id', :fields => ['histology_abbreviation'], :type => :flat, :merge => true, :persist => true, :persist_dir => PCAWG::PROJECT_VAR_DIR.cache
    index.chunked_values_at(specimens).flatten
  end

  property :tumor_specimens => :single do
    PCAWG.donor2tumor_specimen[self]
  end

  property :normal_specimens => :single do
    PCAWG.donor2normal_specimen[self]
  end

  property :SNV_sample => :single do
    Sample.setup(PCAWG.donor2SNV_sample[self], :cohort => cohort).extend AnnotatedArray
  end

  property :CNV_sample => :single do
    Sample.setup(PCAWG.donor2CNV_sample[self], :cohort => cohort).extend AnnotatedArray
  end

  property :expression_samples => :array do
    index = PCAWG.donor2expression_samples
    samples = index.chunked_values_at(self).compact
    Sample.setup(samples, :cohort => cohort)
  end

  property :abbr_color => :single do
    abbr = PCAWG.donor_histology(self, 'histology_abbreviation')
    PCAWG.abb_color(abbr)
  end

  property :gene_status => :single2array do |genes|
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
    self.expression_samples.any?
  end

  property :has_genotype? => :single do
    self.SNV_sample.any?
  end

  property :has_cnv? => :single do
    self.CNV_sample.any?
  end
end

