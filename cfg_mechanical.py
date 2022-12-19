'''
Configuration (cfg) script for Spinal Dorsal Horn Network Model
from Sekiguchi et al. (2021)
'''

from netpyne import specs
from neuron import h
import genrn

cfg = specs.SimConfig()  

cfg.hParams = {'celsius': 36, 'v_init': -60 }
cfg.vrest = cfg.hParams['v_init']
cfg.duration = 5000
 
cfg.recordStep = 0.025

#*--------------------------------*#
#*--- PARAMETERS FOR NETPARAMS ---*#
#*--------------------------------*#

### STIMULATION RATIO OF INPUT FIBERS ###
cfg.stim_ratios = 0.25 # 50mN -> 0.125, 100mN -> 0.25, 200mN -> 0.5, 400mN -> 1.0
if cfg.stim_ratios == 1.0: cfg.freq = '400mN'
elif cfg.stim_ratios == 0.5: cfg.freq = '200mN'
elif cfg.stim_ratios == 0.25:cfg.freq = '100mN'
else: cfg.freq = '50mN'

# SYNAPTIC WEIGHT
cfg.Ab_EX_AMPA = 0.0221559
cfg.Ab_EX_NMDA = 0.015
cfg.Ab_IN_AMPA = 0.00208312
cfg.Ab_IN_NMDA = 0.0098189
cfg.Ad_AMPA = 1.50e-05                
cfg.Ad_NMDA = 1.50e-05                
cfg.C_EX_AMPA = 0.00067115
cfg.C_EX_NMDA = 0.000407345
cfg.C_TrC_AMPA = 0.164705             
cfg.C_DYN_AMPA = 0.184135             
cfg.C_ISLET_AMPA = 0.123535             
cfg.C_NK1_AMPA = 0.00009
cfg.C_NK1_NMDA = 8.7447e-05
cfg.C_NK1_NK1 = 3.2414e-08

cfg.VGLUT3_PKC_AMPA = 0.16629
cfg.VGLUT3_PKC_NMDA = 0.15549
cfg.PV_GABA = 0.29416             
cfg.PV_GLY =  0.011521           
cfg.DYN_ISLET_GABA = 0.36182     
cfg.ISLET_GABA = 0.34293         
cfg.DYN_EX_GABA = 4.50e-05        
cfg.DYN_EX_GLY = 4.50e-05         
cfg.PKC_AMPA = 0.0021
cfg.PKC_NMDA = 0.00315
cfg.TrC_AMPA = 0.00225             
cfg.TrC_NMDA = 0.003              
cfg.VGLUT3_SOM_AMPA = 0.00006
cfg.VGLUT3_SOM_NMDA = 0.00006
cfg.DOR_AMPA = 0.002250          
cfg.DOR_NMDA = 0.002250           
cfg.EX_NK1_AMPA = 8.82981e-06          
cfg.EX_NK1_NMDA = 2.6699e-05
cfg.EX_NK1_NK1 = 9.2715e-07
cfg.DYN_NK1_GABA = 6.3720e-06
cfg.DYN_NK1_GLY = 2.3608e-06

cfg.recordTraces['vs'] = {'sec':'soma', 'loc':0.5,'var':'v'}

# SAVING
cfg.simLabel = '100mN'
cfg.saveFolder = 'data'
cfg.savePickle = False
cfg.saveJson = True
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams']

# ANALYSIS AND PLOTTING
cells = [x for x in range(400, 410, 1)]
cfg.analysis['plotRaster'] = {'include': ['all'], 'timeRange': [0, cfg.duration],'orderInverse': True, 'saveFig': True, 'showFig': True} 
cfg.analysis['plotConn'] = {'includePre': ['all'], 'includePost': ['all'], 'feature': 'weight','logPlot': True, 'saveFig': True, 'showFig': True}
# cfg.analysis['plotSpikeHist'] = {'include': ['eachPop'], 'timeRange': [0,cfg.duration], 'spikeHistBin': 5, 'saveFig': True, 'showFig': False}
# cfg.analysis['plotSpikeStats'] = {'include': ['eachPop'], 'timeRange': [0,cfg.duration], 'saveFig': True, 'showFig': False}
# cfg.analysis['plotTraces'] = {'include': cells, 'timeRange': [0, cfg.duration], 'saveFig': True, 'showFig': False}
# cfg.analysis['plot2Dnet'] = False 

# use for GA simulation
cfg.verbose = True
cfg.filename = 'sim_'
cfg.printPopAvgRates = [ 0, 5000 ]
cfg.dt = 0.025