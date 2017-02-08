module PCAWG

  #PCAWG.claim PCAWG.pandrugs, :proc do |real|
  #  file = DATA_DIR['pandrugs.tar.gz'].find
  #  TmpFile.with_file do |directory|
  #    FileUtils.mkdir_p directory
  #    Misc.in_dir directory do
  #      CMD.cmd("tar xvfz '#{file}'")
  #      Path.setup(directory)

  #      sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
  #      FileUtils.mkdir_p real.find unless File.exists? real.find
  #      directory.glob("pandrugs/*/*.tsv").each do |file|
  #        cohort, donor, *rest = File.basename(file).split("_")
  #        FileUtils.cp file, real[donor].find
  #      end
  #    end
  #    FileUtils.rm_rf directory
  #  end
  #  nil
  #end

  PCAWG.claim PCAWG.pandrugs, :proc do |real|
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        header = nil
        current = nil
        stream = nil
        seen = Set.new
        file = DATA_DIR['pandrugs.tsv.gz'].find
        TSV.traverse file, :type => :array, :bar => "Processing Pandrugs" do |line|
          if line =~ /^#/
            header = "#" << line.split("\t")[2..-1] * "\t"
            next
          end

          cohort, donor, *parts = line.split("\t")

          if current != donor
            stream.close if stream
            stream = Open.open(donor, :mode => "a")
            if not seen.include? donor
              seen << donor
              stream.puts header
            end
          end

          stream.puts parts * "\t"
        end
        stream.close if stream
      end
      FileUtils.mv directory, real
      nil
    end
  end

  PCAWG.claim PCAWG.pandrugs2abbr, :proc do 
    DATA_DIR["tumor_type-Pandrugs-PCAWG.tsv"].tsv :sep2 => /--NONE--/, :type => :single
  end

  PCAWG.claim PCAWG.abbr2pandrugs, :proc do 
    tsv = DATA_DIR["tumor_type-Pandrugs-PCAWG.tsv"].tsv :sep2 => /, */, :type => :flat, :merge => true, :key_field => "tumor_PCAWG"
    ppp tsv.to_s
    tsv.to_s
  end
end
