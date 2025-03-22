import nest
nest.Install('ctcmodule')

# This is the number of readout neurons. Needed in several places below
num_o = 3

nest.resolution = 1
nest.flexible_data_event_size = num_o   # Event size here must be exactly number of readout neurons

# Create neurons:
# - ctc_readout_neuron replaces eprop_readout_neuron_bsshslm_2020
# - additional ctc_loss neuron (only one)
r = nest.Create('eprop_iaf_bsshslm_2020', 10);
o = nest.Create('ctc_readout_neuron', num_o)
t = nest.Create('step_rate_generator')
l = nest.Create('ctc_loss')
                   
# Standard connections for classification task
nest.Connect(r, r, syn_spec='eprop_synapse_bsshslm_2020')
nest.Connect(r, o, syn_spec='eprop_synapse_bsshslm_2020')
nest.Connect(o, r, syn_spec='eprop_learning_signal_connection_bsshslm_2020')
nest.Connect(o, o, syn_spec={'synapse_model': 'rate_connection_delayed', 'receptor_type': 1})
nest.Connect(t, o, syn_spec={'synapse_model': 'rate_connection_delayed', 'receptor_type': 2})

# Connections to and from ctc_loss neuron
# - Connections readout -> ctc_loss
#     - must be rate_connection_instantaneous
#     - must have receptor_types in strictly increasing order given as [[1, 2, ...]]
# - Connections ctc_loss -> readout
#     - must be flexible_data_connection
#     - must be made with same ctc_loss and readout node collections as connection up,
#       just source and target roles swapped
#     - must have same receptor_types as for upward connection, but as [[1], [2], ...]
rports_up = [range(1, len(o)+1)]
rports_down = list(zip(*rports_up))

nest.Connect(o, l, syn_spec={'synapse_model': 'rate_connection_instantaneous',
                                 'receptor_type': rports_up})
nest.Connect(l, o, syn_spec={'synapse_model': 'flexible_data_connection', 'receptor_type': rports_down})

nest.Simulate(4)
