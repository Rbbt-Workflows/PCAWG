require File.join(File.expand_path(File.dirname(__FILE__)),  'test_helper.rb')
require 'pcawg'
require 'rbbt/entity/pcawg'

Log.severity = 0
class TestPCAWG < Minitest::Test
  def _test_genotyped_samples
    study = Study.setup("Bladder-TCC")
    assert study.genotyped_samples.any?
    assert_equal 23, study.genotyped_samples.length
    assert_equal study.genotyped_samples.length, study.genotyped_samples.specimen.length
    assert_equal study.genotyped_samples.specimen.length, study.genotyped_samples.specimen.donor.length
    assert_equal study.genotyped_samples.length, study.genotyped_samples.donor.length
  end

  def test_genotyped_samples_bone
    study = Study.setup("Bone-Osteosarc")
    assert study.genotyped_samples.any?
    assert_equal 46, study.genotyped_samples.length
    assert_equal study.genotyped_samples.length, study.genotyped_samples.specimen.length
    assert_equal study.genotyped_samples.specimen.length, study.genotyped_samples.specimen.donor.length
    assert_equal study.genotyped_samples.length, study.genotyped_samples.donor.length
  end

  def _test_donor_samples
    study = Study.setup("Lymph-CLL")
    assert_equal study.genotyped_samples.first, study.genotyped_samples.first.donor.SNV_sample.first
  end

  def _test_speciment_histology
    study = Study.setup("Bladder-TCC")
    assert_equal ["ENDODERM"], study.genotyped_samples.specimen.histology("histology_tier1").uniq
  end
end

