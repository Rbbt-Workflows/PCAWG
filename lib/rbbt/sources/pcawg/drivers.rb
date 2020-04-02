
module PCAWG

  DRIVER_CALL_FIELD="Brown_observed"

  def self.driver_dir(study)
    clean_study = study.downcase.gsub(/[^a-z]+/, '.').gsub('.tumors','')
    PCAWG.drivers.glob("*").select do |dir| 
      clean_dir = File.basename(dir).downcase.gsub(/[^a-z]+/, '.').gsub('pancan', 'pancancer')
      clean_dir == clean_study
    end.first
  end

  def self.candidate_driver_dir(study)
    clean_study = study.downcase.gsub(/[^a-z]+/, '.').gsub('.tumors','')
    PCAWG.candidate_drivers.glob("*").select do |dir| 
      clean_dir = File.basename(dir).downcase.gsub(/[^a-z]+/, '.').gsub('pancan', 'pancancer')
      clean_dir == clean_study
    end.first
  end

  PCAWG.claim DATA_DIR['final_integration_results_2017_03_16.tar.gz'], :url, "https://b2drop.eudat.eu/s/PWj3NmRJg6KMkX4/download"

  PCAWG.claim PCAWG.drivers, :proc do |real|
    file = DATA_DIR['final_integration_results_2017_03_16.tar.gz'].produce.find
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

  PCAWG.claim PCAWG.donor_drivers_old, :proc do |real|
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

  PCAWG.claim PCAWG.donor_drivers, :proc do |real|
    file = DATA_DIR['pcawg_whitelist_somatic_driver_mutations_beta.csv.gz'].find

    sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD

    fields = []
    fields << "Genomic Mutation"
    fields << "Type"
    fields << "Associated Gene Name"
    fields << "Driver mutation catergory"
    fields << "Driver element catergory"
    fields << "VAF"
    dumper = TSV::Dumper.new :key_field => "icgc_donor_id", :fields => fields, :type => :double, :namespace => PCAWG.organism
    dumper.init
    TSV.traverse file, :header_hash => "", :type => :list, :into => dumper, :fix => Proc.new{|l| l.gsub("|", ";") } do |id, values|

      target, gene, ttype, sample, chr, pos, ref, alt, mut_type, driver_mut_category, driver_element_category, vaf = values
      pos, muts = Misc.correct_vcf_mutation pos.to_i, ref, alt
      res = []
      donor = sample2donor[sample]
      muts.each do |mut|
        mutation = [chr, pos.to_s, mut] * ":"

        res << [donor, [mutation, target, gene, driver_mut_category, driver_element_category, vaf]]
      end
      res.extend MultipleResult
      res
    end
  end

  PCAWG.claim PCAWG.candidate_drivers, :proc do |real|
    excel = DATA_DIR['PCAWG-2,5,9,14 candidate driver release 4-22-2017.xlsx'].find

    require 'rbbt/tsv/excel'
    %w(CDS Promoter 5UTR 3UTR enhancer lncrna.ncrna lncrnaPromoter mirna smallRNA).each_with_index do |name, num|
      tsv = TSV.excel2tsv(excel, :sheet => num, :merge => true)

      tissue_entries = {}
      TSV.traverse tsv do |id, values|
        target, tissues, gene, pvalues, qvalues, num_mutations, num_patients, comments = values
        genes = gene.collect{|g| g.split(";")}.flatten.uniq

        tissues.each_with_index do |tissue, i|
          pvalue = pvalues[i]
          qvalue = qvalues[i]
          num_mutation = num_mutations[i].to_i
          num_patient = num_patients[i].to_i
          comment = (comments || [])[i]

          tissue_entries[tissue] ||= []
          tissue_entries[tissue] << [genes, pvalue, qvalue, num_mutation, num_patient, comment, id]
        end
      end
      
      fields = []
      fields << "P-value"
      fields << "Q-value"
      fields << "Num Mutations"
      fields << "Num Patients"
      fields << "Comment"
      fields << "ID"
      tissue_entries.each do |tissue, values|
        tissue_dir = real[tissue].find
        FileUtils.mkdir_p tissue_dir unless File.exists? tissue_dir

        tissue_values = TSV.setup({}, :key_field => "Associated Gene Name", :fields => fields, :type => :list, :namespace => PCAWG.organism)
        values.each do |list|
          genes, *rest = list
          genes.each do |gene|
            tissue_values[gene] = rest
          end
        end

        Open.write(tissue_dir[name], tissue_values.to_s)
      end
    end
  end

end
