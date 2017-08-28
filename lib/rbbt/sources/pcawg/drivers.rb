
module PCAWG

  DRIVER_CALL_FIELD="Brown_observed"

  def self.driver_dir(study)
    clean_study = study.downcase.gsub(/[^a-z]+/, '.').gsub('.tumors','')
    PCAWG.drivers.glob("*").select do |dir| 
      clean_dir = File.basename(dir).downcase.gsub(/[^a-z]+/, '.').gsub('pancan', 'pancancer')
      clean_dir == clean_study
    end.first
  end

  PCAWG.claim PCAWG.drivers, :proc do |real|
    file = DATA_DIR['final_integration_results_2017_03_16.tar.gz'].find
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        directory.glob("**/*.automatic_method_removal.txt").each do |file|
          cohort, type, *rest = File.basename(file).split(".")
          dumper = TSV::Dumper.new :key_field => "Element", :fields => ["p-value"], :type => :single, :namespace => PCAWG.organism
          dumper.init
          TSV.traverse Open.open(file), :type => :single, :header_hash => "", :into => dumper, :fields => [DRIVER_CALL_FIELD] do |elem, pvalue|
            pvalue = nil if pvalue == "NA"
            [elem, [pvalue]]
          end

          FileUtils.mkdir_p real[cohort].find unless real[cohort].find.exists?
          Misc.consume_stream(dumper.stream, false, real[cohort][type].find)
        end
      end
      FileUtils.rm_rf directory
    end
    nil
  end

  PCAWG.claim PCAWG.drivers_AD, :proc do |real|
    file = Rbbt.root.data.final['ActiveDriver2_AllScores_240117.tgz'].find
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        directory.glob("*/*.observed.txt").each do |file|
          cohort, type, wc, *rest = File.basename(file).split(".")
          type += "." + wc unless wc =~ /Active/
          dumper = TSV::Dumper.new :key_field => "Element", :fields => ["p-value"], :type => :single, :namespace => PCAWG.organism
          dumper.init
          TSV.traverse Open.open(file), :type => :single, :header_hash => "", :into => dumper do |elem, pvalue|
            [elem, [pvalue]]
          end

          FileUtils.mkdir_p real[cohort].find unless real[cohort].find.exists?
          Misc.consume_stream(dumper.stream, false, real[cohort][type].find)
        end
      end
      FileUtils.rm_rf directory
    end
    nil
  end

  PCAWG.claim PCAWG.donor_drivers, :proc do |real|
    file = DATA_DIR['pcawg_whitelist_coding_drivers_v1_sep302016.txt.gz'].find
    
    sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD

    fields = []
    fields << "Genomic Mutation"
    fields << "Associated Gene Name"
    fields << "Protein change"
    fields << "Driver gene"
    fields << "Validated"
    fields << "Source"
    dumper = TSV::Dumper.new :key_field => "icgc_donor_id", :fields => fields, :type => :double, :namespace => PCAWG.organism
    dumper.init
    TSV.traverse file, :header_hash => "", :type => :list, :into => dumper do |sample, values|
      tumor_id, chr, pos, ref, alt, gene, protein, protein_change, consequence, driver_mutation_status, driver_mutation_statement, driver_gene, driver_gene_status, driver_gene_source = values
      pos, muts = Misc.correct_vcf_mutation pos.to_i, ref, alt
      res = []
      donor = sample2donor[sample]
      muts.each do |mut|
        mutation = [chr, pos.to_s, mut] * ":"

        res << [donor, [mutation, gene, protein, driver_gene, driver_gene_status, driver_gene_source]]
      end
      res.extend MultipleResult
      res
    end
  end

end
