from lsst.rapid.analysis.rubinTv import MountTorquePlotter
from lsst.rapid.analysis.utils import checkRubinTvExternalPackages

checkRubinTvExternalPackages()
print('Running mount torque plotter...')
repoDir = '/project/shared/auxTel/'
mountTorquePlotter = MountTorquePlotter(repoDir)
mountTorquePlotter.run()
