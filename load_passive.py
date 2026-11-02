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
data_file = 'D:/infraslow-dynamics/04_data_analysis/004_ibl_bwm/bwmPreprocessedData.pkl'
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
  subject_id = ''
  eid = eids[iExp]

  # Find passive recording period
  try:
    passive_times = one.load_dataset(eid, '*passivePeriods*', collection='alf')
  except:
    pprint('Skipping experiment ' + str(iExp))
    continue
  SP_times = passive_times['spontaneousActivity']
  if len(SP_times):
    
    # Probe 0
    probe_name = 'probe00'
    datasets = one.list_datasets(eid)
    if f'alf/{probe_name}/pykilosort/spikes.times.npy' in datasets:
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
    
    # Probe 1
    probe_name = 'probe01'
    if f'alf/{probe_name}/pykilosort/spikes.times.npy' in datasets:
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
    left_camera = dict()
    if 'alf/_ibl_leftCamera.features.pqt' in datasets:
      video_features = one.load_object(eid, 'leftCamera', collection='alf')
      if hasattr(video_features, 'times'):
        frame_mask = np.logical_and(video_features.times>=SP_times[0],
                                    video_features.times<=SP_times[1])
        if len(frame_mask) <= len(video_features.features.pupilDiameter_raw):
          left_camera = dict(pupilDiameter_raw=video_features.features.pupilDiameter_raw[frame_mask],
                            pupilDiameter_smooth=video_features.features.pupilDiameter_smooth[frame_mask],
                            times=video_features.times[frame_mask])

    # Load right eye pupil diameter
    right_camera = dict()
    if 'alf/_ibl_rightCamera.features.pqt' in datasets:
      video_features = one.load_object(eid, 'rightCamera', collection='alf')
      if hasattr(video_features, 'times'):
        frame_mask = np.logical_and(video_features.times>=SP_times[0],
                                    video_features.times<=SP_times[1])
        if len(frame_mask) <= len(video_features.features.pupilDiameter_raw):
          right_camera = dict(pupilDiameter_raw=video_features.features.pupilDiameter_raw[frame_mask],
                              pupilDiameter_smooth=video_features.features.pupilDiameter_smooth[frame_mask],
                              times=video_features.times[frame_mask])

    # Load wheel movement data
    wheel_data = dict()
    if 'alf/_ibl_wheel.position.npy' in datasets:
      wheel_position = one.load_object(eid, 'wheel', collection='alf')
      wheel_moves = one.load_object(eid, 'wheelMoves', collection='alf')
      wheel_moves.move_intervals = np.empty([0,2])
      for iInterval in range(len(wheel_moves.intervals)):
        interval = wheel_moves.intervals[iInterval,:]
        if interval[1] > SP_times[0] and interval[0] < SP_times[1]:
          interval = np.array([[max(SP_times[0], interval[0]), min(SP_times[1], interval[1])]])
          wheel_moves.move_intervals = np.append(wheel_moves.move_intervals,
                                                 interval, axis=0)
      if hasattr(wheel_position, 'timestamps'):
        frame_mask = np.logical_and(wheel_position.timestamps>=SP_times[0],
                                    wheel_position.timestamps<=SP_times[1])
        if len(frame_mask) <= len(wheel_position.position):
          wheel_data = dict(position=wheel_position.position[frame_mask],
                            times=wheel_position.timestamps[frame_mask],
                            move_intervals=wheel_moves.move_intervals)

    # Store experiment (session) data in a single experiment container
    experiment_data.append(dict(probe0=probe0,
                                probe1=probe1,
                                left_camera=left_camera,
                                right_camera=right_camera,
                                wheel=wheel_data,
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