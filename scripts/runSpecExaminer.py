from lsst.rapid.analysis.rubinTv import SpecExaminer
from lsst.rapid.analysis.utils import checkRubinTvExternalPackages

checkRubinTvExternalPackages()
print('Running spec examiner...')
repoDir = '/project/shared/auxTel/rerun/quickLook'
specExaminer = SpecExaminer(repoDir)
specExaminer.run()
