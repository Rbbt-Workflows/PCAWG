Workflow.require_workflow "Sample"

module Sample

  extend Entity::Identified
  add_identifiers PCAWG.donor_samples, PCAWG::DONOR_FIELD

  helper :watson do
    true
  end

  helper :organism do
    PCAWG.organism
  end

  task :organism => :string do
    PCAWG.organism
  end

  task :genomic_mutations => :array do
    raise "Genomic Mutations not accessible for #{clean_name} due to access limitations"
  end

  property :donor => :array do
    Donor.setup(PCAWG.sample2donor.chunked_values_at self)
  end

  property :specimen => :array do
    Specimen.setup(PCAWG.sample2specimen.chunked_values_at self)
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

Sample.update_task_properties
