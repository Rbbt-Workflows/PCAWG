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

  PCAWG.claim PCAWG.matrices.copy_number, :proc do 
    file = DATA_DIR['all_data_by_genes.rmcnv.pt_170207.txt.gz'].produce
    sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
    donors = [] 
    pos = []
    TSV.traverse file, :type => :array, :into => :stream, :bar => true do |line|
      if line =~ /^Gene Symb/
        fields = line.split("\t")
        fields[0] = "Associated Gene Name"
        aliqs = fields
        fields.each_with_index do |aliq,i| 
          donor = sample2donor[aliq]; 
          if donor
            donors << donor
            pos << i
          end
        end
        "#: organism=" << PCAWG.organism + "\n" +
        "#" + (["Associated Gene Name"] + donors) * "\t"
      else
        all_values = line.split("\t")
        values = all_values.values_at *pos
        
        gene = all_values.first
        gene = gene.split("|chr").first

        ([gene] + values) * "\t"
      end
    end
  end

  PCAWG.claim PCAWG.matrices.copy_number_focal, :proc do 
    file = DATA_DIR['focal_data_by_genes.rmcnv.pt_170207.txt.gz'].produce
    sample2donor = PCAWG.donor_wgs_samples.index :target => PCAWG::DONOR_FIELD
    donors = [] 
    pos = []
    TSV.traverse file, :type => :array, :into => :stream, :bar => true do |line|
      if line =~ /^Gene Symb/
        fields = line.split("\t")
        fields[0] = "Associated Gene Name"
        aliqs = fields
        fields.each_with_index do |aliq,i| 
          donor = sample2donor[aliq]; 
          if donor
            donors << donor
            pos << i
          end
        end
        "#: organism=" << PCAWG.organism + "\n" +
        "#" + (["Associated Gene Name"] + donors) * "\t"
      else
        all_values = line.split("\t")
        values = all_values.values_at *pos
        
        gene = all_values.first
        gene = gene.split("|chr").first

        ([gene] + values) * "\t"
      end
    end
  end
end
