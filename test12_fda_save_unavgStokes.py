# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 18:55:43 2021

Explore how un-averaged Stokes vectors can be saved for PS processing analysis.

@author: vlab_eye
"""

# Import required system libraries for file management
import sys,importlib,os

pc_path = 'C:\\Users\\vlab_eye\\Documents\\local_repo' 
# pc_path = 'C:\\Users\\UCL-SPARC\\Documents\\GitHub'
# pc_path = 'C:\\Users\\cr-eye\\Documents\\GitHub'
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
path_dir = os.path.join(os.getcwd(),'example_data',
                        '[p.Chicken052621][s.SANsChicken_smallV][05-26-2021_13-11-59]')
b_change_editsetting = True
data = Load(directory = path_dir)
data.loadFringe(frame=2)

# sys.exit()
#%% Tomogram processing : complex tomogram, k-space fringe, stokes vectors  
data.reconstructionSettings['processState'] = 'struct+ps' #'angio+ps+hsv+oa'#'+kspace+stokes'
if b_change_editsetting:
    data.reconstructionSettings['spectralBinning'] = True
    data.reconstructionSettings['depthIndex'] = [0,0] # [1100,1500]#(phantom) 
    data.reconstructionSettings['binFract'] = 3
    # data.reconstructionSettings['demodSet'] = [0.4, 0.0, 1.0, 0.0, 0.0, 0.0]
    # data.reconstructionSettings['zoomFactor']= 4
    
if b_change_editsetting:
    data.processOptions['maskOutput'] = False
    data.processOptions['generateProjections'] = False
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
try:
    print("outtom['sv1'].shape >> ", outtom['sv1'].shape)
    print("outtom['iquv1'].shape >> ", outtom['iquv1'].shape)
except:
    print("Error for sv1 or iquv1")
if outtom['k1'] is not None:
    print("outtom['k1'].shape >> ", outtom['k1'].shape)

if outtom['iquv2'] is not None:
    temp_iquv = np.log10(np.squeeze(outtom['iquv2'][:,:,0,3]))
    fig, ax = plt.subplots(1,1)
    ax.imshow(temp_iquv, cmap='gray')
    ax.set(title='unaveraged ii2')
    plt.show()
#%% Structure processing
if b_change_editsetting:
    data.structureSettings['contrastLowHigh'] = [-50, 160] #[0,195]# 
struct_obj = Structure(mode='log')
struct_out = struct_obj.reconstruct(data=data)
for key,val in struct_out.items():
    data.processedData[key] = struct_out[key]
    
print("struct_out.keys() >> ", struct_out.keys())
plt.imshow(cp.asnumpy(struct_out['struct']), aspect ='auto', cmap='gray')


#%% PS processing
if b_change_editsetting:
    data.psSettings['zOffset'] = 15 # this is deltaZ for differential calculation
    data.psSettings['oopFilter']  = 2
    data.psSettings['xFilter'] = 11
    data.psSettings['zFilter']  = 1
    data.psSettings['dopThresh'] = 0.9
    data.psSettings['maxRet'] = 100
    data.psSettings['binFract'] = data.reconstructionSettings['binFract']
    data.psSettings['thetaOffset'] = 0# -int(151/256*180)
    data.psSettings['keepRawStokes'] = True
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


fig, ax = plt.subplots(1,1)
ax.imshow(data.processedData['null'], cmap='gray', aspect='auto')#, clim=(0,2))
plt.show()

# sys.exit()
#%%
# processor = Post()
# processed = processor.processFrameRange(data, procState='struct+ps+stokes+oa+iquv+null', #procState='struct+ps+hsv+oa',#, procAll=True, writeState=True)
#                             procAll=False, startFrame=2, endFrame=2, 
#                             writeState=True,
#                             hold=True, holdState='struct+dop+ret+oa')
# print("\n--------------------")
# print(f"Process one single frame using {type(processor)}, object processed.")
# print("processed.keys() holds the returned results specified by 'holdState' param.")
# print(processed.keys())
# print("--------------------\n")

# # The startFrame should be equal to endFrame in order to show RGB optic axis. 
# fig = plt.figure(figsize=(12, 12))
# ax = fig.add_subplot(221)
# ax.imshow(processed['dop'], cmap='gray', aspect='auto')
# ax.set(title='dop')
# ax = fig.add_subplot(222)
# ax.imshow(processed['ret'], cmap='gray', aspect='auto')
# ax.set(title='ret')
# ax = fig.add_subplot(223)
# ax.imshow(processed['struct'], cmap='gray', aspect='auto')
# ax.set(title='struct')
# ax = fig.add_subplot(224)
# ax.imshow(np.squeeze(processed['oa']/255), aspect='auto')
# ax.set(title='oa') # processed['oa'] ranges from 0 to 1
# plt.show()

# # sys.exit()
#%% Save settings used by Post processor, that were adjusted thus far on Spyder.
print("\n--------------------")
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
print("--------------------\n")

# sys.exit()
#%% Set up the batch processing
# path_directory_project = r'H:\Fig3'
path_directory_project = os.path.split(path_dir)[0]
print(f"b_change_editsetting={b_change_editsetting}.")
b_print_log = True

#------------------------------------------batch processing or single selection
# for session_name in os.listdir(path_directory_project):
#     path_dir = os.path.join(path_directory_project, session_name)
#     print(path_dir)
#------------------------------------------
if True:
#------------------------------------------batch processing or single selection

    if os.path.isdir(path_dir):
        data = Load(directory = path_dir)
        
        # Ensure that the data object still contains the settings parameters defined above.
        if b_change_editsetting:
            for settings_name in list_settings_names:
                setattr(data, settings_name, temp_dict_copy[settings_name]) 
                
        data.reconstructionMode = {'tom': 'heterodyne', 
                                    'struct': 'log', 
                                    'angio': 'cdv', 
                                    'ps': 'sb'}
    
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
        processor.processFrameRange(data, procState='struct+ps+stokes+oa+iquv+null',
                                    procAll=False, startFrame=1, endFrame=10, 
                                    writeState=True)
    else:
        print("Path to a hidden file not a directory")