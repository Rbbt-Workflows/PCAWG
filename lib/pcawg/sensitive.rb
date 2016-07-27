
require 'rbbt/workflow'

Workflow.require_workflow "Sample"

module Sample
  def self.vcf_files(sample)
    filename = sample.split(":").last + '.vcf.gz'
    [PCAWG::PROJECT_DIR.genotypes[filename].find]
  end

  returns "Genomic Mutation"
  task :genomic_mutations => :array do
    io = TSV.get_stream(Sequence.job(:genomic_mutations, sample, :vcf_file => Sample.vcf_files(sample).first).produce)
    CMD.cmd('egrep ":" | sort -u | sed "s/^M:/MT:/" | env LC_ALL=C sort -k1,1 -k2,2n -t:', :in => io, :pipe => true, :no_fail => true)
  end
end

