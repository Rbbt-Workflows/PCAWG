module PCAWG

  PCAWG.claim PCAWG.pandrugs, :proc do |real|
    file = DATA_DIR['pandrugs.tar.gz'].find
    TmpFile.with_file do |directory|
      FileUtils.mkdir_p directory
      Misc.in_dir directory do
        CMD.cmd("tar xvfz '#{file}'")
        Path.setup(directory)

        sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
        FileUtils.mkdir_p real.find unless File.exists? real.find
        directory.glob("pandrugs/*/*.tsv").each do |file|
          cohort, donor, *rest = File.basename(file).split("_")
          FileUtils.cp file, real[donor].find
        end
      end
      FileUtils.rm_rf directory
    end
    nil
  end
end
