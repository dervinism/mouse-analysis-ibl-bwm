% Commonly used data analysis parameters

% IO
processedDataFolder = fullfile('.', 'dependencies');
figFolder = fullfile('.', 'figures');

% Pupil correlation parameters
downsampledRate = 0.2;

% coherence analyses parameters
fRef = 0.03; % Hz
maxFreq = 40; %120;
maxFreq_ca = 40;
maxFreq_pupil = 2.79396772384644;
maxFreq_motion = 2.79396772384644;
winfactor = 4; %10;
freqfactor = 2; %1.6; %1.333;
tapers = 3; %5;
FOI = [40 30 20 15 10 8 6 5 4 3 2 1 0.7 0.5 0.3 0.2 0.1...
  0.07 0.05 0.03 0.02 0.01]; % frequencies of interest (Hz)

% exclusion radius around a unit when calculating local population firing rate
exclRad = 60; %um

% Recorded brain areas
areaLabels = { ...
  'AA',     'AAA',      'Anterior amygdalar area',                                             'association';
  'nCx',    'ACAd',  	  'Anterior cingulate area, dorsal part',                                'association';
  'nCx',    'ACAv',     'Anterior cingulate area, ventral part',                               'association';
  'BG',     'ACB',      'Nucleus accumbens',                                                   '';
  'Th',     'AD',       'Anterodorsal nucleus',                                                'association';
  'Hyp',    'AHN',      'Anterior hypothalamic nucleus',                                       '';
  'paCx',   'AId',      'Agranular insular area, dorsal part',                                 'association';
  'paCx',   'AIp',      'Agranular insular area, posterior part',                              'association';
  'paCx',   'AIv',      'Agranular insular area, ventral part',                                'association';
  'Th',     'AM',       'Anteromedial nucleus',                                                'association';
  'MBr',    'ANcr1',    'Crus 1',                                                              '';
  'MBr',    'ANcr2',    'Crus 2',                                                              '';
  'pCx',    'AON',      'Anterior olfactory nucleus',                                          'sensory';
  'MBr',    'APN',      'Anterior pretectal nucleus',                                          '';
  'nCx',    'APr',      'Area prostriata',                                                     'association';
  'MBr',    'AT',       'Anterior tegmental nucleus',                                          '';
  'nCx',    'AUDd',     'Dorsal auditory area',                                                'association';
  'nCx',    'AUDp',     'Primary auditory area',                                               'sensory';
  'nCx',    'AUDpo',    'Posterior auditory area',                                             'association';
  'nCx',    'AUDv',     'Ventral auditory area',                                               'association';
  'Th',     'AV',       'Anteroventral nucleus of thalamus',                                   'association';
  'Hyp',    'AVP',      'Anteroventral preoptic nucleus',                                      '';
  'AA',     'BLA',      'Basolateral amygdalar nucleus',                                       '';
  'AA',     'BMA',      'Basomedial amygdalar nucleus',                                        '';
  'AA',     'BST',      'Bed nuclei of the stria terminalis',                                  '';
  'Hp',     'CA1',      'Field CA1',                                                           'association';
  'Hp',     'CA2',      'Field CA2',                                                           'association';
  'Hp',     'CA3',      'Field CA3',                                                           'association';
  'AA',     'CEA',      'Central amygdalar nucleus',                                           '';
  'Cereb',  'CENT2',    'Lobule II',                                                           '';
  'Cereb',  'CENT3',    'Lobule III',                                                          '';
  'Th',     'CL',       'Central lateral nucleus of the thalamus',                             'association';
  'CLA',    'CLA',      'Claustrum',                                                           '';
  'BS',     'CLI',      'Central linear nucleus raphe',                                        '';
  'Th',     'CM',       'Central medial nucleus of the thalamus',                              'association';
  'AA',     'COAa',     'Cortical amygdalar area, anterior part',                              '';
  'AA',     'COAp',     'Cortical amygdalar area, posterior part',                             '';
  'Cereb',  'COPY',     'Copula pyramidis',                                                    '';
  'BG',     'CP',       'Caudoputamen',                                                        '';
  'BS',     'CS',       'Superior central nucleus raphe',                                      '';
  'BS',     'CU',       'Cuneate nucleus',                                                     '';
  'Cereb',  'CUL4 5',   'Lobules IV-V',                                                        '';
  'BS',     'CUN',      'Cuneiform nucleus',                                                   '';
  'BS',     'DCO',      'Dorsal cochlear nucleus',                                             '';
  'Hp',     'DG',       'Dentate gyrus',                                                       'association';
  'Hyp',    'DMH',      'Dorsomedial nucleus of the hypothalamus',                             '';
  'Cereb',  'DN',  	    'Dentate nucleus',                                                     '';
  'nCx',    'DP',       'Dorsal peduncular area',                                              'association';
  'BS',     'DR',       'Dorsal nucleus raphe',                                                '';
  'BS',     'DT',       'Dorsal terminal nucleus of the accessory optic tract',                '';
  'BS',     'DTN',      'Dorsal tegmental nucleus',                                            '';
  'paCx',   'ECT',      'Ectorhinal area',                                                     'association';
  'BS',     'ECU',      'External cuneate nucleus',                                            '';
  'paCx',   'ENTl',     'Entorhinal area, lateral part',                                       'association';
  'paCx',   'ENTm',     'Entorhinal area, medial part, dorsal zone',                           'association';
  'AA',     'EPd',      'Endopiriform nucleus, dorsal part',                                   '';
  'AA',     'EPv',      'Endopiriform nucleus, ventral part',                                  '';
  'Th',     'Eth',      'Ethmoid nucleus of the thalamus',                                     'association';
  'Hp',     'FC',       'Fasciola cinerea',                                                    '';
  'Cereb',  'FL',  	    'Flocculus',                                                           '';
  'Cereb',  'FN',       'Fastigial nucleus',                                                   '';
  'nCx',    'FRP',      'Frontal pole, cerebral cortex',                                       'association';
  'BG',     'FS',       'Fundus of striatum',                                                  '';
  'BG',     'GPe',      'Globus pallidus, external segment',                                   '';
  'BS',     'GRN',      'Gigantocellular reticular nucleus',                                   '';
  'nCx',    'GU',       'Gustatory areas',                                                     'sensory';
  'AA',     'HATA',  	  'Hippocampo-amygdalar transition area',                                '';
  'Th',     'IAD',      'Interanterodorsal nucleus of the thalamus',                           'association';
  'Th',     'IAM',      'Interanteromedial nucleus of the thalamus',                           'association';
  'MBr',    'IC',       'Inferior colliculus',                                                 '';
  'Cereb',  'ICB',      'Infracerebellar nucleus',                                             '';
  'Th',     'IGL',      'Intergeniculate leaflet of the lateral geniculate complex',           'sensory';
  'nCx',    'ILA',      'Infralimbic area',                                                    'association';
  'Th',     'IMD',      'Intermediodorsal nucleus of the thalamus',                            'association';
  'Th',     'IP',       'Interposed nucleus',                                                  'association';
  'Cereb',  'IPN',      'Interpeduncular nucleus',                                             '';
  'BS',     'IRN',      'Intermediate reticular nucleus',                                      '';
  'BS',     'ISN',      'Inferior salivatory nucleus',                                         '';
  'Th',     'IntG',     'Intermediate geniculate nucleus',                                     'association';
  'AA',     'LA',       'Lateral amygdalar nucleus',                                           '';
  'BS',     'LAV',      'Lateral vestibular nucleus',                                          '';
  'BS',     'LC',       'Locus ceruleus',                                                      '';
  'Th',     'LD',       'Lateral dorsal nucleus of thalamus',                                  'association';
  'BS',     'LDT',      'Laterodorsal tegmental nucleus',                                      '';
  'Th',     'LGd',      'Dorsal part of the lateral geniculate complex',                       'sensory';
  'Th',     'LGv',      'Ventral part of the lateral geniculate complex',                      'sensory';
  'EPT',    'LH',       'Lateral habenula',                                                    '';
  'Hyp',    'LHA',      'Lateral hypothalamic area',                                           '';
  'nCx',    'LING',     'Lingula (I)',                                                         'association';
  'Th',     'LP',  	    'Lateral posterior nucleus of the thalamus',                           'association';
  'Hyp',    'LPO',      'Lateral preoptic area',                                               '';
  'BS',     'LRN',      'Lateral reticular nucleus',                                           '';
  'BF',     'LSc',      'Lateral septal nucleus, caudal (caudodorsal) part',                   '';
  'BF',     'LSr',      'Lateral septal nucleus, rostral (rostroventral) part',                '';
  'BF',     'LSv',      'Lateral septal nucleus, ventral part',                                '';
  'BS',     'MA',       'Magnocellular nucleus',                                               '';
  'BS',     'MARN',     'Magnocellular reticular nucleus',                                     '';
  'Th',     'MD', 	    'Mediodorsal nucleus of thalamus',                                     'association';
  'BS',     'MDRN',     'Medullary reticular nucleus',                                         '';
  'AA',     'MEA',      'Medial amygdalar nucleus',                                            '';
  'Th',     'MG',       'Medial geniculate complex',                                           'sensory';
  'EPT',    'MH',       'Medial habenula',                                                     '';
  'Hyp',    'MM',       'Medial mammillary nucleus',                                           '';
  'nCx',    'MOp',      'Primary motor area',                                                  'motor';
  'nCx',    'MOs',      'Secondary motor area',                                                'association';
  'Hyp',    'MPN',      'Medial preoptic nucleus',                                             '';
  'Hyp',    'MPO',      'Medial preoptic area',                                                '';
  'MBr',    'MPT',      'Medial pretectal area',                                               '';
  'MBr',    'MRN',      'Midbrain reticular nucleus',                                          '';
  'BF',     'MS',  	    'Medial septal nucleus',                                               '';
  'BS',     'MV',       'Medial vestibular nucleus',                                           '';
  'MBr',    'NB',       'Nucleus of the brachium of the inferior colliculus',                  '';
  'BF',     'NDB',      'Diagonal band nucleus',                                               '';
  'BS',     'NI',       'Nucleus incertus',                                                    '';
  'BS',     'NLL',      'Nucleus of the lateral lemniscus',                                    '';
  'Cereb',  'NOD',      'Nodulus (X)',                                                         '';
  'MBr',    'NOT',      'Nucleus of the optic tract',                                          '';
  'EPT',    'NPC',      'Nucleus of the posterior commissure',                                 '';
  'BS',     'NTS',      'Nucleus of the solitary tract',                                       '';
  'MBr',    'OP',       'Olivary pretectal nucleus',                                           '';
  'nCx',    'ORBm',     'Orbital area, medial part',                                           'association';
  'nCx',    'ORBvl',    'Orbital area, ventrolateral part',                                    'association';
  'BG',     'OT',       'Olfactory tubercle',                                                  'sensory';
  'BS',     'P5',       'Peritrigeminal zone',                                                 '';
  'AA',     'PA',       'Posterior amygdalar nucleus',                                         '';
  'MBr',    'PAG',      'Periaqueductal gray',                                                 '';
  'paCx',   'PAR',      'Parasubiculum',                                                       'association';
  'BS',     'PARN',     'Parvicellular reticular nucleus',                                     '';
  'BS',     'PAS',      'Parasolitary nucleus',                                                '';
  'BS',     'PB',       'Parabrachial nucleus',                                                '';
  'MBr',    'PBG',      'Parabigeminal nucleus',                                               '';
  'BS',     'PC5',      'Parvicellular motor 5 nucleus',                                       '';
  'BS',     'PCG',      'Pontine central gray',                                                '';
  'Th',     'PCN',      'Paracentral nucleus',                                                 'association';
  'BS',     'PDTg',     'Posterodorsal tegmental nucleus',                                     '';
  'paCx',   'PERI',     'Perirhinal area',                                                     'association';
  'Th',     'PF',       'Parafascicular nucleus',                                              'association';
  'Cereb',  'PFL',      'Paraflocculus',                                                       '';
  'BS',     'PG',       'Pontine gray',                                                        '';
  'BS',     'PGRN',     'Paragigantocellular reticular nucleus',                               '';
  'Hyp',    'PH',       'Posterior hypothalamic nucleus',                                      '';
  'Th',     'PIL',      'Posterior intralaminar thalamic nucleus',                             'association';
  'pCx',    'PIR',      'Piriform area',                                                       'sensory';
  'nCx',    'PL',       'Prelimbic area',                                                      'association';
  'Hyp',    'PMd',      'Dorsal premammillary nucleus',                                        '';
  'Th',     'PO',  	    'Posterior complex of the thalamus',                                   'association';
  'Th',     'POL',      'Posterior limiting nucleus of the thalamus',                          'association';
  'paCx',   'POST',  	  'Postsubiculum',                                                       'association';
  'BS',     'PP',       'Peripeduncular nucleus',                                              '';
  'BS',     'PPN',      'Pedunculopontine nucleus',                                            '';
  'MBr',    'PPT',      'Posterior pretectal nucleus',                                         '';
  'Th',     'PR',       'Perireunensis nucleus',                                               'association';
  'paCx',   'PRE',      'Presubiculum',                                                        'association';
  'Cereb',  'PRM',      'Paramedian lobule',                                                   '';
  'BS',     'PRNc',     'Pontine reticular nucleus, caudal part',                              '';
  'BS',     'PRNr',     'Pontine reticular nucleus',                                           '';
  'BS',     'PRP',      'Nucleus prepositus',                                                  '';
  'Hyp',    'PS',       'Parastrial nucleus',                                                  '';
  'Hyp',    'PSTN',     'Parasubthalamic nucleus',                                             '';
  'BS',     'PSV',      'Principal sensory nucleus of the trigeminal',                         '';
  'Th',     'PT',       'Parataenial nucleus',                                                 'association';
  'Hyp',    'PVH',      'Paraventricular hypothalamic nucleus',                                '';
  'Hyp',    'PVHd',     'Paraventricular hypothalamic nucleus, descending division',           '';
  'Th',     'PVT',      'Paraventricular nucleus of the thalamus',                             'association';
  'Cereb',  'PYR',      'Pyramus (VIII)',                                                      '';
  'BS',     'Pa5',      'Paratrigeminal nucleus',                                              '';
  'Hyp',    'PeF',      'Perifornical nucleus',                                                '';
  'Th',     'PoT',      'Posterior triangular thalamic nucleus',                               'association';
  'paCx',   'ProS',     'Prosubiculum',                                                        'association';
  'Th',     'RE',       'Nucleus of reuniens',                                                 'association';
  'Th',     'RH',       'Rhomboid nucleus',                                                    'association';
  'BS',     'RM',       'Nucleus raphe magnus',                                                '';
  'MBr',    'RN',       'Red nucleus',                                                         '';
  'MBr',    'RPF',      'Retroparafascicular nucleus',                                         '';
  'MBr',    'RR',       'Midbrain reticular nucleus, retrorubral area',                        '';
  'paCx',   'RSPagl',   'Retrosplenial area, lateral agranular part',                          'association';
  'paCx',   'RSPd',     'Retrosplenial area, dorsal part',                                     'association';
  'paCx',   'RSPv',     'Retrosplenial area, ventral part',                                    'association';
  'Th',     'RT',       'Reticular nucleus of the thalamus',                                   'association';
  'MBr',    'SCm',      'Superior colliculus, motor related',                                  '';
  'MBr',    'SCs',      'Superior colliculus, sensory related',                                '';
  'BF',     'SF',       'Septofimbrial nucleus',                                               '';
  'Th',     'SGN',      'Suprageniculate nucleus',                                             'association';
  'BF',     'SI',       'Substantia innominata',                                               '';
  'Cereb',  'SIM',      'Simple lobule',                                                       '';
  'Th',     'SMT',      'Submedial nucleus of the thalamus',                                   'association';
  'MBr',    'SNc',      'Substantia nigra, compact part',                                      '';
  'MBr',    'SNr',      'Substantia nigra, reticular part',                                    '';
  'BS',     'SOC',      'Superior olivary complex',                                            '';
  'Th',     'SPA',      'Subparafascicular area',                                              'association';
  'Th',     'SPF',      'Subparafascicular nucleus',                                           'association';
  'BS',     'SPIV',     'Spinal vestibular nucleus',                                           '';
  'BS',     'SPVC',     'Spinal nucleus of the trigeminal, caudal part',                       '';
  'BS',     'SPVI',     'Spinal nucleus of the trigeminal, interpolar part',                   '';
  'BS',     'SPVO',     'Spinal nucleus of the trigeminal, oral part',                         '';
  'nCx',    'SSp-bfd',  'Primary somatosensory area, barrel field',                            'sensory';
  'nCx',    'SSp-ll',   'Primary somatosensory area, lower limb',                              'sensory';
  'nCx',    'SSp-m',    'Primary somatosensory area, mouth',                                   'sensory';
  'nCx',    'SSp-n',    'Primary somatosensory area, nose',                                    'sensory';
  'nCx',    'SSp-tr',   'Primary somatosensory area, trunk',                                   'sensory';
  'nCx',    'SSp-ul',   'Primary somatosensory area, upper limb',                              'sensory';
  'nCx',    'SSp-un',   'Primary somatosensory area, unassigned',                              'sensory';
  'nCx',    'SSs',      'Supplemental somatosensory area',                                     'association';
  'STh',    'STN',      'Subthalamic nucleus',                                                 '';
  'paCx',   'SUB',      'Subiculum',                                                           'association';
  'Hyp',    'SUM',      'Supramammillary nucleus',                                             '';
  'BS',     'SUT',      'Supratrigeminal nucleus',                                             '';
  'BS',     'SUV',      'Superior vestibular nucleus',                                         '';
  'Th',     'SubG',     'Subgeniculate nucleus',                                               'association';
  'nCx',    'TEa',      'Temporal association areas',                                          'association';
  'BS',     'TRN',      'Tegmental reticular nucleus',                                         '';
  'BF',     'TRS',      'Triangular nucleus of septum',                                        '';
  'pCx',    'TTd',      'Taenia tecta, dorsal part',                                           'sensory';
  'pCx',    'TTv',      'Taenia tecta, ventral part',                                          'sensory';
  'Hyp',    'TU',       'Tuberal nucleus',                                                     '';
  'Cereb',  'UVU',      'Uvula (IX)',                                                          '';
  'BS',     'V',        'Motor nucleus of trigeminal',                                         '';
  'Th',     'VAL',      'Ventral anterior-lateral complex of the thalamus',                    'motor';
  'BS',     'VCO',      'Ventral cochlear nucleus',                                            '';
  'BS',     'VII',      'Facial motor nucleus',                                                '';
  'nCx',    'VISC',     'Visceral area',                                                       'sensory';
  'nCx',    'VISa',     'Anterior area',                                                       'association';
  'nCx',    'VISam',    'Anteromedial visual area',                                            'association';
  'nCx',    'VISl',     'Lateral visual area',                                                 'association';
  'nCx',    'VISli',    'Laterointermediate area',                                             'association';
  'nCx',    'VISp',     'Primary visual area',                                                 'sensory';
  'nCx',    'VISpl',    'Posterolateral visual area',                                          'association';
  'nCx',    'VISpm',    'posteromedial visual area',                                           'association';
  'paCx',   'VISpor',   'Postrhinal area',                                                     'association';
  'nCx',    'VISrl',    'Rostrolateral visual area',                                           'association';
  'Th',     'VM',       'Ventral medial nucleus of the thalamus',                              'association';
  'Hyp',    'VMH',      'Ventromedial hypothalamic nucleus',                                   '';
  'Th',     'VPL',      'Ventral posterolateral nucleus of the thalamus',                      'sensory';
  'Th',     'VPLpc',    'Ventral posterolateral nucleus of the thalamus, parvicellular part',  'sensory';
  'Th',     'VPM',      'Ventral posteromedial nucleus of the thalamus',                       'sensory';
  'Th',     'VPMpc',    'Ventral posteromedial nucleus of the thalamus, parvicellular part',   'sensory';
  'MBr',    'VTA',      'Ventral tegmental area',                                              '';
  'MBr',    'VTN',      'Ventral tegmental nucleus',                                           '';
  'Cereb',  'VeCB',     'Vestibulocerebellar nucleus',                                         '';
  'Th',     'Xi',       'Xiphoid thalamic nucleus',                                            'association';
  'STh',    'ZI',       'Zona incerta',                                                        '';
  '???',    'root',     'root',                                                                '';
  '???',    'void',     'void',                                                                '';
  'BS',     'x',        'Nucleus x',                                                           ''};

% Selected areas of interest
areasOI = {'sensory'; 'association'; 'motor'; 'Th'; 'nCx'; 'paCx'; 'pCx'; ...
  'Hp'; 'sensory-Th'; 'association-Th'; 'motor-Th'; 'sensory-nCx'; ...
  'association-nCx'; 'motor-nCx'; 'association-paCx'; 'sensory-pCx'; ...
  'VPL-VPM'; 'LG'; 'PO'; 'LP'; 'VPL-VPM-LG-PO-LP'; 'SSp-bfd'; 'SSp-body'; ...
  'SSp'; 'RSPagl'; 'RSPd'; 'RSPv'; 'RSPd-RSPv'; 'RSP'; 'CA1'; 'CA2'; 'CA3'; ...
  'DG'; 'CA-DG'};

% Minimal area set
% areasMinimal = {'Th'; 'nCx'; 'VPL-VPM-LG-PO-LP'; 'SSp-bfd'; 'SSp-body'; ...
%   'SSp'; 'RSPagl'; 'RSPd'; 'RSPv'; 'RSPd-RSPv'; 'RSP'; 'CA1'; 'CA2'; ...
%   'CA3'; 'CA'; 'DG'; 'CA-DG'};
areasMinimal = {'VPL-VPM'; 'VPL-VPM-LG'; 'VPL-VPM-LG-PO'; 'VPL-VPM-LG-PO-LP'; ...
  'SSp-bfd'; 'SSp'; 'RSPagl'; 'RSP'; 'CA'; 'DG'; 'CA-DG'};