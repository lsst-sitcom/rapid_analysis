from lsst.rapid.analysis.rubinTv import ImExaminer
from lsst.rapid.analysis.utils import checkRubinTvExternalPackages

checkRubinTvExternalPackages()
print('Running imExaminer...')
repoDir = '/project/shared/auxTel/rerun/quickLook'
imExaminer = ImExaminer(repoDir)
imExaminer.run()
