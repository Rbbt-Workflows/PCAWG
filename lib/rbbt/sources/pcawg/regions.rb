module PCAWG
  PCAWG.claim PCAWG.enhancer_ranges, :proc do
    file = PCAWG::DATA_DIR.annotations.map["map.enhancer.gene"].find
    TSV.traverse file, :type => :array, :into => :stream do |line|
      range, genes = line.split("\t")
      chr,start,eend = range.split(/:|-/)
      chr.sub!('chr','')
      [chr, start, eend, genes.split(";").uniq*";"] * "\t"
    end
  end

  PCAWG.claim PCAWG.regions, :proc do |directory|
    FileUtils.mkdir_p directory.find
    PCAWG::DATA_DIR.annotations.gene.glob("*.bed").each do |file|
      name = File.basename file
      io = Misc.sort_mutation_stream(Open.open(file), "\t")
      Misc.consume_stream(io, false, directory[name].find)
    end
  end
end
