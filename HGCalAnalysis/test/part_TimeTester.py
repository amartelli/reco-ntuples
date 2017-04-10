import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


import FWCore.Utilities.FileUtils as FileUtils
readFiles = cms.untracked.vstring()
#readFiles.extend(FileUtils.loadListFromFile ('/afs/cern.ch/work/a/amartell/HGCAL/clustering/CMSSW_9_0_0_pre2/src/HGCalCalibration/HitValidation/test/fileList/listKL130_Pt100.txt') )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        "/store/cmst3/group/hgcal/CMG_studies/Production/partGun_amartell_PDGid22_nPart1_E60_900pre2_20170113/RECO/partGun_PDGid22_x50_E60.0To60.0_RECO_99.root" 
        ),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )


from RecoNtuples.HGCalAnalysis.timeRecHitEstimator_cfi import HGCalTimeEstimator

process.ana = cms.EDAnalyzer('HGCTimeTester',
                             HGCEEInput = cms.InputTag('HGCalRecHit:HGCEERecHits'),
                             HGCFHInput = cms.InputTag('HGCalRecHit:HGCHEFRecHits'),
                             HGCBHInput = cms.InputTag('HGCalRecHit:HGCHEBRecHits'),
                             dEdXweights = HGCalTimeEstimator.dEdXweights,
                             thicknessCorrection = HGCalTimeEstimator.thicknessCorrection,
                             HGCEE_fCPerMIP = HGCalTimeEstimator.HGCEE_fCPerMIP,
                             HGCEE_noisefC = HGCalTimeEstimator.HGCEE_noisefC,
                             HGCEF_noisefC = HGCalTimeEstimator.HGCEF_noisefC,
                             HGCBH_noiseMIP = HGCalTimeEstimator.HGCBH_noiseMIP,
                             HGCEE_keV2fC  = HGCalTimeEstimator.HGCEE_keV2fC,
                             HGCHEF_keV2fC = HGCalTimeEstimator.HGCHEF_keV2fC,
                             HGCHB_keV2MIP = HGCalTimeEstimator.HGCHB_keV2MIP
                             )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:HGC_recHitTimeTester.root")
                                   )


process.p = cms.Path(process.ana)

