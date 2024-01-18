module PCAWG

  Rbbt.claim DATA_DIR['October_2016_whitelist_2583.snv_mnv_indel.maf.gz'].tap{|o| o.resource = PCAWG}, :proc do |filename|
    raise "You do not have permission to view Genomic Mutations in this server. Otherwise, please place the file October_2016_whitelist_2583.snv_mnv_indel.maf.gz into #{ filename }"
  end

  PCAWG.claim PCAWG.genotype_counts, :proc do |real|

    TmpFile.with_file do |directory|
      Path.setup(directory)
      aliquote2donor = PCAWG.donor_wgs_samples.tsv :key_field => 'tumor_wgs_aliquot_id', :type => :single

      if Open.exists? directory
        FileUtils.rm_rf directory.find
      end
      FileUtils.mkdir_p directory 
      last_donor = nil
      io = nil

      TSV.traverse DATA_DIR['October_2016_whitelist_2583.snv_mnv_indel.maf.gz'].produce, :type => :array, :bar => "Processing mutations" do |line|
        next if line =~ /^Tumor_Sample_Barcode/
        parts = line.split("\t")
        chr, start, eend, ref, alt, alt2, ali, alt_count, ref_count  = parts.values_at 1,2,3,7,8,9,12,37,38

        donor = aliquote2donor[ali]
        next if donor.nil?
        start = start.to_i
        eend = eend.to_i

        _muts = ([alt, alt2].compact.uniq - [ref])
        _mut  = _muts.first

        pos, muts = Misc.correct_vcf_mutation start, ref, _mut
        mutations = muts.collect{|m| [chr, pos, m] * ":"}

        mutation_vafs = mutations.collect{|m| [m, alt_count, ref_count] * "\t" }

        if last_donor != donor
          io.close if io
          if File.exist?(directory[donor])
            io = Open.open(directory[donor], :mode => 'a')
          else
            io = Open.open(directory[donor], :mode => 'w')
            io.puts "#: :type=:list#:namespace=#{PCAWG.organism}"
            io.puts ["#Genomic Mutation", "Alt. reads", "Ref. reads"] * "\t"
          end
          io.puts(mutation_vafs * "\n")
        else
          io.puts(mutation_vafs * "\n")
        end
      end
      io.close
      FileUtils.mv directory.find, real.find
    end
    nil
  end

  PCAWG.claim PCAWG.genotypes, :proc do |real|
    PCAWG.genotype_counts.produce.find.glob("*").each do |counts_file|
      donor = File.basename(counts_file)
      io = CMD.cmd("cut -f 1 '#{counts_file}'|grep -v '#'")
      Open.write(real[donor], io)
    end
  end
end
