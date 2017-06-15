module PCAWG
  PCAWG.claim PCAWG.clonality.structure, :proc do |real|
    file = DATA_DIR['20170325_consensus_subclonal_reconstruction_beta1.tar.gz'].find
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        directory.glob("*/*_subclonal_structure.txt.gz").each do |file|
          sample = File.basename(file).split("_").first
          donor = sample2donor[sample]
          next if donor.nil?
          tsv = file.tsv :header_hash => "", :type => :list, :cast => :to_f

          Open.write(real[donor].find, tsv.to_s)
        end
      end
      FileUtils.rm_rf directory
    end
    nil
  end

  PCAWG.claim PCAWG.clonality.assignments, :proc do |real|
    file = DATA_DIR['20170325_consensus_subclonal_reconstruction_beta1.tar.gz'].find
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        directory.glob("*/*_cluster_assignments.txt.gz").each do |file|
          sample = File.basename(file).split("_").first
          donor = sample2donor[sample]
          next if donor.nil?

          final = TSV.setup({}, :key_field => "Genomic Position", :fields =>  ["Assignments"], :type => :flat)
          TSV.traverse file, :type => :array, :into => final do |line|
            next if line =~ /chromosome/
            chr, pos, type, *rest = line.split("\t")
            chr2 = rest.pop
            pos2 = rest.pop
            if chr2 == "NA"
              cpos = [chr, pos] * ":"
            else
              cpos = [chr, pos] * ":" + ":-:" + [chr2, pos2] * ":"
            end
            [cpos, rest]
          end

          Open.write(real[donor].find, final.to_s)
        end
      end
      FileUtils.rm_rf directory
    end
    nil
  end

  PCAWG.claim PCAWG.clonality.timing, :proc do |real|
    file = DATA_DIR['20170325_consensus_subclonal_reconstruction_beta1.tar.gz'].find
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        directory.glob("*/*_mutation_timing.txt.gz").each do |file|
          sample = File.basename(file).split("_").first
          donor = sample2donor[sample]
          next if donor.nil?

          final = TSV.setup({}, :key_field => "Genomic Position", :fields =>  ["Timing"], :type => :single)
          TSV.traverse file, :type => :array, :into => final do |line|
            next if line =~ /chromosome/
            chr, pos, type, *rest = line.split("\t")
            pos2 = rest.pop
            chr2 = rest.pop
            if chr2 == "NA"
              cpos = [chr, pos] * ":"
            else
              cpos = [chr, pos] * ":" + ":-:" + [chr2, pos2] * ":"
            end
            [cpos, rest]
          end

          Open.write(real[donor].find, final.to_s)
        end
      end
      FileUtils.rm_rf directory
    end
    nil
  end
end
