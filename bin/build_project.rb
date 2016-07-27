#!/usr/bin/env ruby

require 'rbbt-util'
require 'rbbt/util/simpleopt'

$0 = "rbbt #{$previous_commands*""} #{ File.basename(__FILE__) }" if $previous_commands

options = SOPT.setup <<EOF

Build PCAWG project with vcf files

$ #{$0} [options] <filename.tsv|->

-h--help Print this help

EOF
if options[:help]
  if defined? rbbt_usage
    rbbt_usage 
  else
    puts SOPT.usage
  end
  exit 0
end


release_file = 'data/release_mar2016.v1.tsv'
snv_dir = 'data/snv_mnv'
indel_dir = 'data/indel'
matrix_file = 'data/tophat_star_fpkm_uq.tsv'
project_dir = 'project/'

contaminated_samples = Open.read('data/contaminated_samples').split("\n")
oxo_tsv = TSV.open(release_file, :header_hash => '', :key_field => 'tumor_wgs_aliquot_id', :fields => ['tumor_wgs_oxog_score'], :type => :single, :cast => :to_f, :sep2 => ',')
tsv = TSV.open(release_file, :header_hash => '', :key_field => 'tumor_wgs_aliquot_id', :fields => ['dcc_project_code'], :type => :double, :sep2 => ',').to_single
rnaseq_tsv = TSV.open(release_file, :header_hash => '', :key_field => 'tumor_wgs_aliquot_id', :fields => ['tumor_rna_seq_tophat_alignment_bam_file_name'], :type => :double, :sep2 => ',').to_single
rnaseq_normal_tsv = TSV.open(release_file, :header_hash => '', :key_field => 'tumor_wgs_aliquot_id', :fields => ['normal_rna_seq_tophat_alignment_bam_file_name'], :type => :double, :sep2 => ',').to_single

#- Remove oxyg <= 40
#- 19 contaminated
#- LOWSUPPORT

metadata =<<EOF
---
:organism: Hsa/feb2014
:watson: true
:exome: false
:users:
  - pcawg
:expression_type: fpkm
:condition: [CONDITION]
EOF

matrix_samples = {}
missing_oxo = 0
bad_oxo = 0
contaminated = 0
TSV.traverse tsv, :bar => "Placing vcf" do |id,study|
  if contaminated_samples.include? id
    Log.warn "Sample #{ id } skipped for contamination problems"
    contaminated += 1
    next
  end

  if oxo_tsv[id].nil?
    Log.warn "Sample #{ id } skipped for MISSING oxo score"
    missing_oxo +=1
    next
  end
  if oxo_tsv[id] <= 40
    Log.warn "Sample #{ id } skipped for oxo problems"
    bad_oxo +=1
    next
  end

  study_dir = File.join(project_dir, study)
  genotype_dir = File.join(study_dir, 'genotypes/vcf')
  snv_file = File.join(snv_dir, id + '.merged.somatic.snv_mnv.vcf.gz')
  indel_file = File.join(indel_dir, id + '.merged.somatic.indel.vcf.gz')

  if not File.directory? study_dir
    FileUtils.mkdir_p genotype_dir unless File.directory? genotype_dir
    Open.write(File.join(study_dir, 'metadata.yaml'), metadata.sub('[CONDITION]', study.sub(/-.*/,'')))
  end

  if File.exists? snv_file
    target_file = File.join(genotype_dir, id + '.vcf')
    CMD.cmd("zcat #{ snv_file } | grep -v 'LOWSUPPORT' > #{target_file}", :no_fail => true)
    CMD.cmd("zcat #{ indel_file } | grep -v '^#' | grep -v 'LOWSUPPORT' >> #{target_file}", :no_fail => true)
    CMD.cmd("gzip #{ target_file }")
  end

  if rnaseq_tsv[id] and not rnaseq_tsv[id].empty?
    matrix_samples[study] ||= {}
    matrix_samples[study][:tumor] ||= {}
    matrix_samples[study][:tumor][id] ||= []
    matrix_samples[study][:tumor][id] << rnaseq_tsv[id].split(".")[1]
  end

  if rnaseq_normal_tsv[id] and not rnaseq_normal_tsv[id].empty?
    matrix_samples[study] ||= {}
    matrix_samples[study][:normal] ||= {}
    matrix_samples[study][:normal][id] ||= []
    matrix_samples[study][:normal][id] << rnaseq_normal_tsv[id].split(".")[1]
  end
end

ppp "Missing oxo info: #{ missing_oxo }"
ppp "Bad oxo info: #{ bad_oxo }"
ppp "Contaminated: #{ contaminated }"

all_matrix_samples = TSV.parse_header(matrix_file, :header_hash => '').fields

TSV.traverse matrix_samples, :bar => "Placing rnaseq" do |study,sample_info|
  if sample_info[:tumor]
    study_dir = File.join(project_dir, study)
    study_matrix_file = File.join(study_dir, 'matrices/tumor_rnaseq/data')
    FileUtils.mkdir_p File.dirname(study_matrix_file) unless File.exists? File.dirname(study_matrix_file)

    tumor_samples = []
    sample_names = []
    sample_info[:tumor].each do |id,mids|
      mid = mids.first
      next unless all_matrix_samples.include? mid
      tumor_samples << mids.first
      sample_names << id
    end
    tumor_sample_positions = tumor_samples.collect{|s| all_matrix_samples.index(s)+1}

    CMD.cmd("cut -f 1,#{tumor_sample_positions*","} #{matrix_file} | sed 's/^feature/#Ensembl Gene ID/;s/^\\(ENSG[[:digit:]]*\\)\.[[:digit:]]*/\\1/' > #{study_matrix_file}")
  end

  if sample_info[:normal]
    study_dir = File.join(project_dir, study)
    study_matrix_file = File.join(study_dir, 'matrices/normal_rnaseq/data')
    FileUtils.mkdir_p File.dirname(study_matrix_file) unless File.exists? File.dirname(study_matrix_file)

    tumor_samples = []
    sample_names = []
    sample_info[:normal].each do |id,mids|
      mid = mids.first
      next unless all_matrix_samples.include? mid
      tumor_samples << mids.first
      sample_names << id
    end
    tumor_sample_positions = tumor_samples.collect{|s| all_matrix_samples.index(s)+1}

    CMD.cmd("cut -f 1,#{tumor_sample_positions*","} #{matrix_file} | sed 's/^feature/#Ensembl Gene ID/;s/^\\(ENSG[[:digit:]]*\\)\.[[:digit:]]*/\\1/' > #{study_matrix_file}")
  end
end
