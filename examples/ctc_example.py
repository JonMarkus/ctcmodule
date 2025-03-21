import nest
nest.Install('ctcmodule')
nest.resolution = 1
nest.wfr_interpolation_order = 3

r = nest.Create('eprop_iaf_bsshslm_2020', 10);
o = nest.Create('ctc_readout_neuron', 3)
t = nest.Create('step_rate_generator')
l = nest.Create('ctc_loss')

# standard connections for classification task
nest.Connect(r, r, syn_spec='eprop_synapse_bsshslm_2020')
nest.Connect(r, o, syn_spec='eprop_synapse_bsshslm_2020')
nest.Connect(o, r, syn_spec='eprop_learning_signal_connection_bsshslm_2020')
nest.Connect(o, o, syn_spec={'synapse_model': 'rate_connection_delayed', 'receptor_type': 1})
nest.Connect(t, o, syn_spec={'synapse_model': 'rate_connection_delayed', 'receptor_type': 2})

# ctc connection
rports_up = [range(1, len(o)+1)]
rports_down = list(zip(*rports_up))

nest.Connect(o, l, syn_spec={'synapse_model': 'rate_connection_instantaneous',
                                 'receptor_type': rports_up})
#nest.Connect(l, o, syn_spec={'synapse_model': 'gap_junction', 'receptor_type': rports_down})

nest.Simulate(10)
