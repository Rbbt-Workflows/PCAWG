module PCAWG

  def self.concept_color(concept)
    @@color_palette ||= PCAWG.color_palette.tsv
    key = concept.gsub(/[^a-zA-Z]/,'.').downcase
    @@color_palette[key]
  end

  def self.abbr_color(abbr)
    concept_color(abbr)
  end

  #def self.abbr_color(abbr)
  #  @@abbr_color ||= PCAWG.subtype_info.tsv :fields => ["Color (RGB code)"], :type => :single
  #  @@abbr_color[abbr]
  #end

  #def self.abbr_colors
  #  PCAWG.subtype_info.tsv :fields => ["Color (RGB code)"], :type => :single
  #end

  def self.abbr_driver_genes
    @@abbr_driver_genes ||= PCAWG.driver_genes.tsv :fields => ["Ensembl Gene ID"], :type => :flat
  end

  def self.all_driver_genes
    @@all_driver_genes ||= PCAWG.driver_genes.tsv(:key_field => "Ensembl Gene ID", :fields => []).keys
  end

  def self.activating_driver_genes
    @@activating_driver_genes ||= PCAWG.driver_genes.tsv(:key_field => "Ensembl Gene ID", :fields => ["Role"]).select("Role" => "Act").keys
  end


  def self.lof_driver_genes
    @@lof_driver_genes ||= PCAWG.driver_genes.tsv(:key_field => "Ensembl Gene ID", :fields => ["Role"]).select("Role" => "LoF").keys
  end

  #{{{ Color palette

  PCAWG.claim PCAWG.color_palette, :proc do
    source_url = "https://raw.githubusercontent.com/ICGC-TCGA-PanCancer/pcawg-colour-palette/master/pcawg.colour.palette.R"
    txt = Open.read(source_url, :nocache => true)
    palette = TSV.setup({}, :key_field => "Concept key", :fields => ["Color"], :type => :single)
    txt.scan(/^\s+([a-z0-9\.]+)\s*<-\s*'(#[a-zA-Z0-9]+)'/m).each do |values|
      key, value = values
      palette[key] = value
    end
    palette.to_s
  end
end
