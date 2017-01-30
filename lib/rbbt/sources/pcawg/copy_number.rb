module PCAWG
  PCAWG.claim PCAWG.CNV, :proc do |real|
    file = DATA_DIR['consensus.20170119.somatic.cna.annotated.tar.gz'].find
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        directory.glob("*.txt").each do |file|
          sample = File.basename(file).split(".").first
          donor = sample2donor[sample]
          next if donor.nil?
          stream = TSV.traverse Open.open(file), :type => :array, :into => :stream do |line|
            next if line.include? 'chrom'
            chr1, start, eend, total, minor, major, stars = line.split("\t")
            next if stars == "NA"
            [chr1, start, eend, total, minor, major, stars] * ":"
          end

          Open.write(real[donor].find, stream.read)
        end
      end
      FileUtils.rm_rf directory
    end
    nil
  end
end
