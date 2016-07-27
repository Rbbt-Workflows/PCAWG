require 'rbbt-util'
require 'rbbt/workflow'
require 'pcawg'

require 'pcawg/sensitive' if false

Misc.add_libdir #if __FILE__ == $0

#require 'rbbt/sources/PCAWG'

module PCAWG
  extend Workflow

end

#require 'PCAWG/tasks/basic.rb'

#require 'rbbt/knowledge_base/PCAWG'
#require 'rbbt/entity/PCAWG'

