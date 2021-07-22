from lsst.rapid.analysis.rubinTv import Monitor
from lsst.rapid.analysis.utils import checkRubinTvExternalPackages

checkRubinTvExternalPackages()
print('Running monitor...')
repoDir = '/project/shared/auxTel/rerun/quickLook'
monitor = Monitor(repoDir)
monitor.run()
