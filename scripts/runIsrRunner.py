from lsst.rapid.analysis.rubinTv import IsrRunner

print('Running isr runner...')
repoDir = '/project/shared/auxTel/'
isrRunner = IsrRunner(repoDir)
isrRunner.run()
