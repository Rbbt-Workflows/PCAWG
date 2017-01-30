module PCAWG
  PCAWG.claim PCAWG.SV, :proc do |real|
    file = DATA_DIR['pcawg_consensus_1.6.161116.somatic_svs.tar.gz'].find
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        directory.glob("*.bedpe.gz").each do |file|
          sample = File.basename(file).split(".").first
          donor = sample2donor[sample]
          next if donor.nil?
          stream = TSV.traverse Open.open(file), :type => :array, :into => :stream do |line|
            next if line.include? 'chrom'
            chr1, pos1, _pos1, chr2, pos2, _pos2, _id, support, strand1, strand2, sv_class = line.split("\t")
            [chr1, pos1, sv_class,chr2,pos2] * ":"
          end

          Open.write(real[donor].find, stream.read)
        end
      end
      FileUtils.rm_rf directory
    end
    nil
  end
end
