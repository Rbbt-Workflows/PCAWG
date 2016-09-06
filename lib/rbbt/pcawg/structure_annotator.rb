module PCAWG
  def self.PCAWG_residues
    @@PCAWG_residues ||= Persist.persist_tsv(nil, "PCAWG::residues", {}, :persist => true, :serializer => :list, :dir => PCAWG::PROJECT_VAR_DIR.Structure) do |data|
                           isoform_residue_mutations = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Mutated Isoform"], :type => :flat)

                           job = Study.setup("PCAWG").mi(:job)
                           job.produce

                           TSV.traverse job, :type => :array do |mi|
                             next unless mi =~ /(ENSP\d+):([A-Z])(\d+)([A-Z])$/ and $2 != $4
                             residue = $3.to_i
                             protein = $1
                             isoform_residue_mutations[[protein, residue] * ":"] ||= []
                             isoform_residue_mutations[[protein, residue] * ":"] << mi
                           end

                           data.merge! isoform_residue_mutations
                           isoform_residue_mutations.annotate data

                           data
                         end
  end
  def self.PCAWG_mi_annotations
    @@PCAWG_mi_annotations ||= begin
                                 fields = [
                                   PCAWG::DONOR_FIELD,
                                   PCAWG::SNV_SAMPLE_FIELD,
                                   'histology_abbreviation',
                                   'histology_tier1',
                                   'histology_tier2',
                                   'histology_tier3',
                                   'histology_tier4',
                                   'dcc_project_code',
                                 ]
                                 Persist.persist_tsv(nil, "PCAWG::annotations", { :key_field => "Genomic Mutations", :fields => fields}, {:serializer => :double, :persist => true, :dir => PCAWG::PROJECT_VAR_DIR.Structure} ) do |data|
                                   organism = PCAWG.organism
                                   mi_annotations = TSV.setup({}, :key_field => "Mutated Isoform", :fields => fields, :type => :double, :namespace => organism)

                                   job = Study.setup("PCAWG").mi_incidence(:job)
                                   job.produce

                                   sample2donor = PCAWG.donor_samples.index :target => PCAWG::DONOR_FIELD
                                   donor2histology = PCAWG.specimen_histology.tsv :key_field => PCAWG::DONOR_FIELD, :fields => fields.select{|f| f=~/histo/}
                                   donor2code = PCAWG.donor_sample_info.tsv :key_field => PCAWG::DONOR_FIELD, :fields => ["dcc_project_code"], :type => :single

                                   TSV.traverse job, :into => mi_annotations, :bar => "Processing PCAWG mi annotations" do |mi,samples|
                                     mi = mi.first if Array === mi
                                     annotations = []
                                     samples.flatten.each do |sample|
                                       donor = sample2donor[sample]
                                       histologies = donor2histology[donor].collect{|l| l*"&"}
                                       code = donor2code[donor]
                                       annotations << [donor,sample] + histologies + [code]
                                     end
                                     [mi, Misc.zip_fields(annotations)]
                                   end

                                   data.merge! mi_annotations
                                   mi_annotations.annotate data
                                   data
                                 end
                               end
  end
end

Workflow.require_workflow "Structure"

module Structure

    ANNOTATORS["PCAWG"] = Annotator.new PCAWG::DONOR_FIELD,
      PCAWG::SNV_SAMPLE_FIELD,
      'histology_abbreviation',
      'histology_tier1',
      'histology_tier2',
      'histology_tier3',
      'histology_tier4',
      'dcc_project_code' do |isoform, residue,organism|

        @PCAWG_residues ||= PCAWG.PCAWG_residues
        @PCAWG_mi_annotations ||= PCAWG.PCAWG_mi_annotations

        isoform_residue = isoform + ":" << residue.to_s
        mis = @PCAWG_residues[isoform_residue]

        next if mis.nil?
        next if mis.empty?
        mis.uniq!
        tmp = {}
        annots = mis.collect{|mi| @PCAWG_mi_annotations[mi]}
        Misc.zip_fields(annots)
      end
end
