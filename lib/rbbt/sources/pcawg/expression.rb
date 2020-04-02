module PCAWG

  PCAWG.claim DATA_DIR['joint_fpkm_uq.tsv.gz'], :url, "https://b2drop.eudat.eu/s/Ym4WNiwMJd3STgj/download"

  PCAWG.claim PCAWG.matrices.gene_expression, :proc do 
    TSV.traverse DATA_DIR["joint_fpkm_uq.tsv.gz"].produce, :type => :array, :into => :stream, :bar => true do |line|
      if line =~ /^feature/
        fields = line.split("\t")
        fields[0] = "Ensembl Gene ID"
        "#" + fields * "\t"
      else
        values = line.split("\t")
          values[0] = values[0].split(".").first
          values * "\t"
      end
    end
  end
end
