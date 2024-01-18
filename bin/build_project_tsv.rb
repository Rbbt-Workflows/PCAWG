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


snv_file = 'data/tsv/snv.tsv.m2c.sam.black.oxo40.noRECA'
indel_file = 'data/tsv/indel.tsv.m2c.sam.black.oxo40.noRECA'
release_file = 'data/release_mar2016.v1.tsv'

rnaseq_tsv = TSV.open(release_file, :header_hash => '', :key_field => 'tumor_wgs_aliquot_id', :fields => ['tumor_rna_seq_tophat_alignment_bam_file_name'], :type => :double, :sep2 => ',').to_single
rnaseq_normal_tsv = TSV.open(release_file, :header_hash => '', :key_field => 'tumor_wgs_aliquot_id', :fields => ['normal_rna_seq_tophat_alignment_bam_file_name'], :type => :double, :sep2 => ',').to_single

matrix_file = 'data/tophat_star_fpkm_uq.tsv'
project_dir = 'project/'

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

Log.severity = 0
matrix_samples = {}
sample_streams = {}
TSV.traverse snv_file, :bar => "Extracting variants", :type => :array do |line|
  next if line =~ /^CHROM/
  chrom, pos, rsid, ref, alt, _qual, _filter, info, id, dcc_project_code, tumor_wgs_oxog_score, exclude_sample = line.split("\t")
  study = dcc_project_code

  pos, alts = Misc.correct_vcf_mutation(pos,ref,alt)
  muts = alts.collect{|alt| [chrom, pos, alt] * ":" }

  study_dir = File.join(project_dir, study)
  genotype_dir = File.join(study_dir, 'genotypes/vcf')

  if not File.directory? study_dir
    FileUtils.mkdir_p genotype_dir unless File.directory? genotype_dir
    Open.write(File.join(study_dir, 'metadata.yaml'), metadata.sub('[CONDITION]', study.sub(/-.*/,'')))
  end

  sample_streams[id] ||= begin
                         sample_genotype = File.join(genotype_dir, id)
                         Open.open(sample_genotype, :mode => 'w')
                       end

  sample_streams[id].puts muts * "\n"


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

TSV.traverse indel_file, :bar => "Extracting variants", :type => :array do |line|
  next if line =~ /^CHROM/
  chrom, pos, rsid, ref, alt, _qual, _filter, info, id, dcc_project_code, tumor_wgs_oxog_score, exclude_sample = line.split("\t")
  study = dcc_project_code

  pos, alts = Misc.correct_vcf_mutation(pos.to_i,ref,alt)
  muts = alts.collect{|alt| [chrom, pos, alt] * ":" }

  sample_streams[id] ||= begin
                           sample_genotype = File.join(project_dir, study, 'genotypes', id)
                           FileUtils.mkdir_p File.join(project_dir, study, 'genotypes') unless File.exist?(File.join(project_dir, study, 'genotypes'))
                           Open.open(sample_genotype, :mode => 'w')
                         end

  sample_streams[id].puts muts * "\n"
end

all_matrix_samples = TSV.parse_header(matrix_file, :header_hash => '').fields

TSV.traverse matrix_samples, :bar => "Placing rnaseq" do |study,sample_info|
  if sample_info[:tumor]
    study_dir = File.join(project_dir, study)
    study_matrix_file = File.join(study_dir, 'matrices/tumor_rnaseq/data')
    FileUtils.mkdir_p File.dirname(study_matrix_file) unless File.exist? File.dirname(study_matrix_file)

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
    FileUtils.mkdir_p File.dirname(study_matrix_file) unless File.exist? File.dirname(study_matrix_file)

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
