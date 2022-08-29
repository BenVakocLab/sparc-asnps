# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 18:55:43 2021

Run processing for beginning of the FDA data for JNE manuscript
@author: vlab_eye
"""

# Import required system libraries for file management
import sys,importlib,os

pc_path = 'C:\\Users\\vlab_eye\\Documents\\local_repo' #
# pc_path = 'C:\\Users\\UCL-SPARC\\Documents\\GitHub'
# Provide path to oct-cbort library
module_path = os.path.abspath(os.path.join(pc_path, 'oct-cbort'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
# Import oct-cbort library
from oct import *

# ASN libs
module_path2=os.path.abspath(os.path.join(pc_path, 'asnlibs'))
if module_path2 not in sys.path:
    sys.path.append(module_path2)
from mgh_io import MGHio
from asndev_analysis import mark_line
import asnlib_mplset
# Import additional libraries
import shutil, time, copy


#%% Load data
path_dir = r'D:\OFDIData\user.Ilyas\[p.S3B_11_25_19][s._S_loc3][11-25-2019_12-38-31]'
# path_dir = r'H:\Fig4\[p.S4A_12_3_19][s._baseline][12-03-2019_10-15-18]'
# path_dir = r'H:\C1\[p.C1B_11_12_19][s._C_loc2][11-12-2019_13-21-43]'
b_change_editsetting = True
data = Load(directory = path_dir)
data.loadFringe(frame=1300)


#%% Tomogram processing : complex tomogram, k-space fringe, stokes vectors  
data.reconstructionSettings['processState'] = 'struct+angio+ps+hsv+oa'#'+kspace+stokes'
if b_change_editsetting:
    data.reconstructionSettings['spectralBinning'] = True
    data.reconstructionSettings['depthIndex'] = [1000,2000] # [1100,1500]#(phantom) 
    data.reconstructionSettings['binFract'] = 3
    data.reconstructionSettings['demodSet'] = [0.4, 0.0, 1.0, 0.0, 0.0, 0.0]
    data.reconstructionSettings['zoomFactor']= 4
    
if b_change_editsetting:
    data.processOptions['maskOutput'] = False
    data.processOptions['generateProjections'] = True
    data.processOptions['projState'] = 'struct+hsv'
    data.processOptions['projType'] = 'sum'
    data.processOptions['OOPAveraging'] = True
    data.processOptions['correctSystemOA'] = False
    

tom = Tomogram(mode='heterodyne')
outtom = tom.reconstruct(data=data)
for key,val in outtom.items():
    data.processedData[key] = outtom[key]
    
print("outtom.keys() >> ", outtom.keys())
plt.imshow(cp.asnumpy(cp.log10(cp.abs(outtom['tomch1'][:,:]))), aspect ='auto', cmap='gray')
plt.imshow(cp.asnumpy(cp.log10(cp.abs(outtom['tomch2'][:,:]))), aspect ='auto', cmap='gray')

print("outtom['tomch1'].shape >> ", outtom['tomch1'].shape)
print("outtom['sv1'].shape >> ", outtom['sv1'].shape)
if outtom['k1'] is not None:
    print("outtom['k1'].shape >> ", outtom['k1'].shape)

sys.exit()
#%% Structure processing
if b_change_editsetting:
    data.structureSettings['contrastLowHigh'] = [0,195]# [-50, 160]
struct_obj = Structure(mode='log')
struct_out = struct_obj.reconstruct(data=data)
for key,val in struct_out.items():
    data.processedData[key] = struct_out[key]
    
print("struct_out.keys() >> ", struct_out.keys())
plt.imshow(cp.asnumpy(struct_out['struct']), aspect ='auto', cmap='gray')


#%% PS processing
if b_change_editsetting:
    data.psSettings['zOffset'] = 5 # this is deltaZ for differential calculation
    data.psSettings['oopFilter']  = 2
    data.psSettings['xFilter'] = 11
    data.psSettings['zFilter']  = 1
    data.psSettings['dopThresh'] = 0.9
    data.psSettings['maxRet'] = 100
    data.psSettings['binFract'] = data.reconstructionSettings['binFract']
    data.psSettings['thetaOffset'] = 0# -int(151/256*180)
    data.hsvSettings['maskThresholds'] = np.array([230, 15, 45])
    data.hsvSettings['structWeight'] = np.array([40, 70])
    data.hsvSettings['dopWeight'] = np.array([200, 250])
    data.hsvSettings['retWeight'] = np.array([5, 128])
    data.hsvSettings['thetaRef'] = 0#-int(151/256*180)
    

print(data.psSettings)
ps = Polarization('sym') # Polarization('sb')
outps = ps.reconstruct(data=data)
for key,val in outps.items():
    data.processedData[key] = outps[key]
    
print("outps.keys() >> ", outps.keys())

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(221)
ax.imshow(data.processedData['dop'], cmap='gray', aspect='auto')
ax.set(title='dop')
ax = fig.add_subplot(222)
ax.imshow(data.processedData['ret'], cmap='gray', aspect='auto')
ax.set(title='ret')
ax = fig.add_subplot(223)
ax.imshow(data.processedData['theta'], aspect='auto')
ax.set(title='theta')
ax = fig.add_subplot(224)
ax.imshow(data.processedData['oa'], aspect='auto')
ax.set(title='oa')
plt.show()


SVF1, SVF2, SVN1, SVN2, QUV1, QUV2 = filtNormStokes(data.processedData['sv1'],
                                                    data.processedData['sv2'],
                                                    stokesFilter=ps.filter)
#%%
# sys.exit()
processor = Post()
processed = processor.processFrameRange(data, procState='struct+ps+hsv+oa',#, procAll=True, writeState=True)
                            procAll=False, startFrame=1000, endFrame=1001, writeState=False,
                            hold=True, holdState='struct+oa+ret+dop')

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(221)
ax.imshow(processed['dop'], cmap='gray', aspect='auto')
ax.set(title='dop')
ax = fig.add_subplot(222)
ax.imshow(processed['ret'], cmap='gray', aspect='auto')
ax.set(title='ret')
ax = fig.add_subplot(223)
ax.imshow(processed['struct'], cmap='gray', aspect='auto')
ax.set(title='struct')
ax = fig.add_subplot(224)
ax.imshow(processed['oa'], aspect='auto')
ax.set(title='oa')
plt.show()
#%%
if b_change_editsetting:
    print("Settings changed in Spyder")
else:
    print("Saved edit settings used")

# Make temporary dict to save the data settings
list_settings_names = ['reconstructionSettings', 'structureSettings', 'angioSettings', 'psSettings', 'hsvSettings', 'processOptions']
temp_dict_copy = dict.fromkeys(list_settings_names)
for settings_name in list_settings_names:
    temp_dict_copy[settings_name] = copy.deepcopy(getattr(data, settings_name))
    print(f"{settings_name} saved in the temp_dict.")
    
# sys.exit()
#%% Set up the batch processing
# path_directory_project = r'H:\Fig3'
path_directory_project = os.path.split(path_dir)[0]
b_change_editsetting = True
b_print_log = True

#---------------------------------------------
# for session_name in os.listdir(path_directory_project):
#     path_dir = os.path.join(path_directory_project, session_name)
#     print(path_dir)
#---------------------------------------------
if True:
    path_dir = r'D:\OFDIData\user.Ilyas\[p.S3B_11_25_19][s._S_loc3][11-25-2019_12-38-31]'
#---------------------------------------------
    if os.path.isdir(path_dir):
        data = Load(directory = path_dir)
        
        # Ensure that the data object still contains the settings parameters defined above.
        if b_change_editsetting:
            for settings_name in list_settings_names:
                setattr(data, settings_name, temp_dict_copy[settings_name]) 
                
        data.reconstructionMode = {'tom': 'heterodyne', 
                                    'struct': 'log', 
                                    'angio': 'cdv', 
                                    'ps': 'sym'}
    
        # Set up logging to print on IDE console
        if b_print_log:
            logging.basicConfig(stream=sys.stdout, level=logging.INFO)
            logger = logging.getLogger()
            logger.setLevel(logging.INFO)
            logger.__dict__
    
        #%% Process and save the whole volume
        # Write editsettings.ini file and copy to the settings folder
        data.generateEditSettings()
        src = os.path.join(path_dir, 'editsettings.ini')
        dst = os.path.join(path_dir, 'Processed', 'spyder_used_editsettings.ini')
        shutil.copy(src, dst)
    
        # Initialize the post processor. 
        processor = Post()
        processor.processFrameRange(data, procState='struct+ps+stokes+oa',#, procAll=True, writeState=True)
                                    procAll=False, startFrame=1101, endFrame=1101, writeState=True)
    else:
        print("Path to a hidden file not a directory")