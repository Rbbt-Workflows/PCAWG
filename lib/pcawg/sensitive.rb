
require 'rbbt/workflow'

Workflow.require_workflow "Sample"

PCAWG::PROJECT_DIR.genotypes.produce
module Sample
  def self.vcf_files(sample)
    sample = sample.split(":").last
    PCAWG::PROJECT_DIR.genotypes[sample].glob('*.vcf.gz')
  end

  returns "Genomic Mutation"
  task :genomic_mutations => :array do
    ios = Sample.vcf_files(sample).collect do |vcf|
      TSV.get_stream(Sequence.job(:genomic_mutations, sample, :vcf_file => vcf).produce)
    end
    io = Misc.intercalate_streams(ios)
    CMD.cmd('egrep ":" | sort -u | sed "s/^M:/MT:/" | env LC_ALL=C sort -k1,1 -k2,2n -t:', :in => io, :pipe => true, :no_fail => true)
  end
end

