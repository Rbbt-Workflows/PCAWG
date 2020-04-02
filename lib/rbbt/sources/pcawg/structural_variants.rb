module PCAWG

  PCAWG.claim DATA_DIR['pcawg_consensus_1.6.161116.somatic_svs.tar.gz'], :url, "https://b2drop.eudat.eu/s/q2APXgY5HRjYjZC/download"

  PCAWG.claim PCAWG.SV, :proc do |real|
    file = DATA_DIR['pcawg_consensus_1.6.161116.somatic_svs.tar.gz'].produce

    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file.find}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        directory.glob("*.bedpe.gz").each do |file|
          begin
          sample = File.basename(file).split(".").first
          donor = sample2donor[sample]

          next if donor.nil?

          stream = TSV.traverse file, :type => :array, :into => :stream do |line|
            next if line.include? 'chrom'
            chr1, pos1, _pos1, chr2, pos2, _pos2, _id, support, strand1, strand2, sv_class = line.split("\t")
            [chr1, pos1, sv_class,chr2,pos2] * ":"
          end

          Misc.sensiblewrite(real[donor].find, stream)

          rescue
            Log.exception $!
            raise $!
          end
        end
      end
      FileUtils.rm_rf directory
    end
    nil
  end
end
