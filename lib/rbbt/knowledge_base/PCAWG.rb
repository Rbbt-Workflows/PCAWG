
require 'rbbt/knowledge_base'
require 'rbbt/sources/organism'

module PCAWG

  class << self 
    attr_accessor :knowledge_base_dir
  end
  self.knowledge_base_dir = Rbbt.var.knowledge_base.PCAWG

  def self.organism
    Organism.default_code("Hsa")
  end

  def self.knowledge_base
    @knowledge_base ||= begin
                          kb = KnowledgeBase.new self.knowledge_base_dir, self.organism

                          kb
                        end
  end
end

PCAWG.knowledge_base.register :gene_principal_isoform_mutations do
  Workflow.require_workflow "Appris"
  require 'rbbt/sources/appris'

  study = Study.setup("PCAWG")

  all_mis = study.mi(:principal => true)

  protein_genes = Organism.transcripts(PCAWG.organism).tsv :persist => true, :key_field => "Ensembl Protein ID", :fields => ["Ensembl Gene ID"], :type => :single

  mi2genes = Misc.process_to_hash(all_mis) do |mis|
    mis.collect do |mi|
      protein, sep, change = mi.partition ":"
      next unless protein =~ /ENSP/
      next unless Appris.principal_isoform_list.include? protein
      protein_genes[protein]
    end
  end

  tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Mutated Isoform"], :type => :flat, :namespace => PCAWG.organism)

  all_mis.each do |mi|
    gene = mi2genes[mi]
    next if gene.nil?
    tsv[gene] ||= []
    tsv[gene] << mi
  end

  tsv
end


PCAWG.claim PCAWG.gene_damage_analysis, :proc do

  tsv = TSV.setup({}, :key_field => "Ensembl Gene ID",
                  :fields => ["Avg. damage score", "Bg. Avg. damage score", "T-test p-value"],
                  :type => :list, :cast => :to_f, :namespace => PCAWG.organism, :unnamed => true)

  Workflow.require_workflow 'DbNSFP'
  require 'rbbt/util/R'
  database = DbNSFP.database
  database.unnamed = true
  damage_fields = database.fields.select{|f| f =~ /rankscore/ }
  db = PCAWG.knowledge_base.get_database(:gene_principal_isoform_mutations)

  R.eval "a=1" # To start Rserver for all cpus
  RbbtSemaphore.with_semaphore 1 do |sem|
    TSV.traverse db, :cpus => 10, :into => tsv, :bar => "Damage analysis using DbNSFP" do |gene, mis|
      mis = mis.flatten.compact
      next if mis.empty?
      protein = mis.first.partition(":").first
      next unless protein =~ /^ENSP/

      dbNSFP_tsv = database.get_prefix(protein).slice(damage_fields)
      dbNSFP_tsv.unnamed = true

      all_damage_scores = dbNSFP_tsv.collect{|k,values| good = values.reject{|v| v == -999}; good.any? ? Misc.mean(good) : nil}.compact
      damage_scores = dbNSFP_tsv.select(:key => mis).collect{|k,values| good = values.reject{|v| v == -999}; good.any? ? Misc.mean(good) : nil}.compact

      if damage_scores.length < 3 or all_damage_scores.uniq.length < 3
        damage_score_pvalue = 1
      else
        RbbtSemaphore.synchronize(sem) do
          damage_score_pvalue = R.eval("t.test(#{R.ruby2R(damage_scores)}, #{R.ruby2R(all_damage_scores)}, 'greater')['p.value']").to_f
        end
      end
      [gene, [Misc.mean(damage_scores), Misc.mean(all_damage_scores), damage_score_pvalue]]
    end
  end

  tsv.to_s
end
