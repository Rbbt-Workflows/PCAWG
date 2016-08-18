

module Specimen
  extend Entity

  annotation :cohort

  property :donor => :array do
    Donor.setup(PCAWG.specimen2donor.chunked_values_at self)
  end

  property :histology => :array do |key|
    index = PCAWG.specimen_histology_index
    column = index.column(key)
    column.chunked_values_at self
  end

  property :abbreviation => :array do 
    histology('histology_abbreviation')
  end

  property :color => :array do 
    abbreviation.collect{|abbr|
      PCAWG.abb_color(abbr)
    }
  end

  property :has_genotype? => :array do
    donor.has_genotype?
  end

  property :has_cnv? => :array do
    donor.has_cnv?
  end

  property :has_expression? => :array do
    donor.has_expression?
  end
end
