import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = 'MCRUN2_74_V7'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )


process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    '/store/mc/RunIISpring15MiniAODv2/GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/AC1B0462-EB6D-E511-8B83-29EC1E0159C7.root',
#'/store/mc/RunIISpring15MiniAODv2/GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/AA80A3A1-EA6D-E511-933C-002590593920.root',
#'/store/mc/RunIISpring15MiniAODv2/GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/B6289F14-EB6D-E511-ABB9-074059FCAAED.root',
#'/store/mc/RunIISpring15MiniAODv2/GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/F271F9B2-EA6D-E511-9A9A-A5D954A458DE.root'
#'/store/group/phys_higgs/resonant_HH/RunII/MicroAOD/HHSignal74X/1_3_0/GluGluToBulkGravitonToHHTo2B2G_M-300_narrow_13TeV-madgraph/HHSignal74X-1_3_0-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/160201_123855/0000/myMicroAODOutputFile_1.root'
#'/store/mc/RunIISpring15MiniAODv2/GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/AA80A3A1-EA6D-E511-933C-002590593920.root'
'/store/mc/RunIISpring15MiniAODv2/GluGluToBulkGravitonToHHTo2B2G_M-300_narrow_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/B0234D29-056F-E511-9F89-02163E016BCE.root'
    )
)

process.demo = cms.EDAnalyzer('GenAnalyzer')

from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
#from RecoJets.JetProducers.jetToolbox_cff import jetToolbox
process.myJetSequence = cms.Sequence()


process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

#jetToolbox( process, 'ak4', 'myJetSequence', 'outTemp'#, JETCorrPayload='None'#, btagDiscriminators='pfCombinedInclusiveSecondaryVertexV2BJetTags'
#             ) #, addPrunedSubjets=True )


process.demo.outFile = cms.untracked.string("aatest.root")
process.demo.photonTag = cms.InputTag('slimmedPhotons')
process.demo.genTag = cms.InputTag('prunedGenParticles')
process.demo.genjetTag = cms.InputTag('slimmedGenJets')
process.demo.jetTag = cms.InputTag('slimmedJets')


process.p = cms.Path(process.myJetSequence*process.demo)
