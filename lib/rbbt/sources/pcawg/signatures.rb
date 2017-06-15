module PCAWG

  PCAWG.claim PCAWG.mutation_signature.profiles, :proc do 
    file = DATA_DIR['PCAWG_signature_patterns_beta2.20170320.txt.gz'].find
    tsv = file.tsv :key_field => "Mutation Subtype", :type => :list, :header_hash => ''
    good_fields = tsv.fields.select{|f| f.include? "Signature" }
    tsv.slice(good_fields).transpose("Signature").to_s
  end

  PCAWG.claim PCAWG.mutation_signature.assignments, :proc do 
    file = DATA_DIR['PCAWG_sub_signatures_in_samples_beta2.20170320.txt.gz'].find
    tsv = file.tsv :key_field => "aliquot_id", :type => :list, :header_hash => ''
    tsv.key_field = "tumor_wgs_aliquot_id"
    good_fields = tsv.fields.select{|f| f.include? "Signature" }
    tsv = tsv.slice(good_fields).change_key "icgc_donor_id", :identifiers => PCAWG.donor_wgs_samples
    tsv.to_s
  end

  PCAWG.claim PCAWG.mutation_signature.probabilities, :proc do 
    file = DATA_DIR['All_samples_PCAWG_probabilities_for_subs_beta2.20170320.txt.gz'].find
    tsv = file.tsv :key_field => "aliquot_id", :type => :list, :header_hash => ''
    tsv.key_field = "tumor_wgs_aliquot_id"
    good_fields = tsv.fields.select{|f| f.include? "Signature" }
    tsv = tsv.slice(good_fields).change_key "icgc_donor_id", :identifiers => PCAWG.donor_wgs_samples
    tsv.to_s
  end
end
