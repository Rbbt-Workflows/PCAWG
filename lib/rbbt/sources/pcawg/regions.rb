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
end
