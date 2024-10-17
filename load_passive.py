from one.api import ONE
from brainbox.population.decode import get_spike_counts_in_bins
from brainbox.io.one import SpikeSortingLoader, SessionLoader
from brainbox.ephys_plots import plot_brain_regions
#from brainbox.behavior.wheel import velocity
from brainbox.task.trials import get_event_aligned_raster, get_psth
from ibllib.atlas import AllenAtlas
from brainwidemap import bwm_query, load_good_units, load_trials_and_mask, bwm_units
#from brainwidemap.imbizo.encoding_functions import get_choice_time_shuffle
import one.alf.io as alfio
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from pprint import pprint
import os.path
import pickle

# Load data file
data_file = 'C:/Users/44079/Work/Leicester/infraslow-dynamics/04_data_analysis/005_ibl_bwm/bwmPreprocessedData.pkl'
if os.path.isfile(data_file):
  with open(data_file, 'rb') as f:
    bwmPreprocessedData = pickle.load(f)
    experiment_data = bwmPreprocessedData['experiment_data']
else:
  experiment_data = list()

# Set up network connection
ONE.setup(base_url='https://openalyx.internationalbrainlab.org', silent=True)
one = ONE(password='international')

# Load brain area ID mapping
ba = AllenAtlas()

# Generate BWM probe table
probeTable_df = bwm_query()

# Generate BWM unit table
#unit_df = bwm_units(one, rt_range=(0.08, 0.2), min_errors=3, min_qc=1., min_units_sessions=(10, 2))
unitTable_df = bwm_units(one)

# List experiments
eids = probeTable_df.eid.unique()

# Load experiment data
for iExp in range(len(eids)):
  pprint('Progress: ' + str(iExp) + '/' + str(len(eids)))
  eid = eids[iExp]

  # Find passive recording period
  try:
    passive_times = one.load_dataset(eid, '*passivePeriods*', collection='alf')
  except:
    pprint('Skipping experiment ' + str(iExp))
    continue
  SP_times = passive_times['spontaneousActivity']
  if len(SP_times):
    
    # Probe 1
    probe_name = 'probe00'
    collections = one.list_collections(eid)
    if f'alf/{probe_name}/pykilosort' in collections:
      probe_mask = np.logical_and(probeTable_df.eid==eid,
                                  probeTable_df.probe_name==probe_name)
      spikes = one.load_object(eid, 'spikes', collection=f'alf/{probe_name}/pykilosort', 
                              attribute=['times', 'clusters'])
      spike_mask = np.logical_and(spikes.times>=SP_times[0],
                                  spikes.times<=SP_times[1])
      clusters = one.load_object(eid, 'clusters', collection=f'alf/{probe_name}/pykilosort', 
                                attribute=['channels', 'waveforms', 'uuids'])
      channels = one.load_object(eid, 'channels', collection=f'alf/{probe_name}/pykilosort', 
                                attribute=['localCoordinates', 'brainLocationIds_ccf_2017', 'labels'])
      probe0 = dict(spike_times=spikes.times[spike_mask],
                    spike_clusters=spikes.clusters[spike_mask],
                    probe_name=probe_name,
                    pid=probeTable_df.pid[probe_mask],
                    cluster_channels=clusters.channels,
                    cluster_waveforms=clusters.waveforms,
                    cluster_uuids=clusters.uuids,
                    channels_localCoordinates=channels.localCoordinates,
                    channels_brainLocationIds_ccf_2017=channels.brainLocationIds_ccf_2017,
                    channels_labels=channels.labels,
                    channels_brain_areas=ba.regions.id2acronym(
                      channels.brainLocationIds_ccf_2017, mapping='Beryl'))
      subject_id = probeTable_df.subject[probe_mask]
    else:
      probe0 = dict()
    
    # Probe 2
    probe_name = 'probe01'
    if f'alf/{probe_name}/pykilosort' in collections:
      probe_mask = np.logical_and(probeTable_df.eid==eid,
                                  probeTable_df.probe_name==probe_name)
      spikes = one.load_object(eid, 'spikes', collection=f'alf/{probe_name}/pykilosort', 
                               attribute=['times', 'clusters'])
      spike_mask = np.logical_and(spikes.times>=SP_times[0],
                                  spikes.times<=SP_times[1])
      clusters = one.load_object(eid, 'clusters', collection=f'alf/{probe_name}/pykilosort', 
                                attribute=['channels', 'waveforms', 'uuids'])
      channels = one.load_object(eid, 'channels', collection=f'alf/{probe_name}/pykilosort', 
                                attribute=['localCoordinates', 'brainLocationIds_ccf_2017', 'labels'])
      probe1 = dict(spike_times=spikes.times[spike_mask],
                    spike_clusters=spikes.clusters[spike_mask],
                    probe_name=probe_name,
                    pid=probeTable_df.pid[probe_mask],
                    cluster_channels=clusters.channels,
                    cluster_waveforms=clusters.waveforms,
                    cluster_uuids=clusters.uuids,
                    channels_localCoordinates=channels.localCoordinates,
                    channels_brainLocationIds_ccf_2017=channels.brainLocationIds_ccf_2017,
                    channels_labels=channels.labels,
                    channels_brain_areas=ba.regions.id2acronym(
                      channels.brainLocationIds_ccf_2017, mapping='Beryl'))
      subject_id = probeTable_df.subject[probe_mask]
    else:
      probe1 = dict()
    
    # Load left eye pupil diameter
    video_features = one.load_object(eid, 'leftCamera', collection='alf')
    frame_mask = np.logical_and(video_features.times>=SP_times[0],
                                video_features.times<=SP_times[1])
    left_camera = dict(pupilDiameter_raw=video_features.features.pupilDiameter_raw[frame_mask],
                       pupilDiameter_smooth=video_features.features.pupilDiameter_smooth[frame_mask],
                       times=video_features.times[frame_mask])

    # Load right eye pupil diameter
    video_features = one.load_object(eid, 'rightCamera', collection='alf')
    frame_mask = np.logical_and(video_features.times>=SP_times[0],
                                video_features.times<=SP_times[1])
    right_camera = dict(pupilDiameter_raw=video_features.features.pupilDiameter_raw[frame_mask],
                        pupilDiameter_smooth=video_features.features.pupilDiameter_smooth[frame_mask],
                        times=video_features.times[frame_mask])

    # Store experiment (session) data in a single experiment container
    experiment_data.append(dict(probe0=probe0,
                                probe1=probe1,
                                left_camera=left_camera,
                                right_camera=right_camera,
                                spontaneous_activity_times=SP_times,
                                eid=eid,
                                subject_id=subject_id,
                                iteration_state=iExp))

    # Store probe and unit tables
    bwmPreprocessedData = dict(experiment_data=experiment_data,
                               probeTable=probeTable_df,
                               unitTable=unitTable_df,
                               iteration_state=iExp)

# Save data
with open(data_file, 'wb') as f:
  pickle.dump(bwmPreprocessedData, f)
pprint('Data loading complete.')