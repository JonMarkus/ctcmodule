import nest
nest.Install('ctcmodule')
nest.resolution = 1

r = nest.Create('eprop_iaf_bsshslm_2020', 10);
nest.Connect(r, r, syn_spec='eprop_synapse_bsshslm_2020')
o = nest.Create('ctc_readout_neuron', 3)
nest.Connect(r, o, syn_spec='eprop_synapse_bsshslm_2020')
nest.Connect(o, r, syn_spec='eprop_learning_signal_connection_bsshslm_2020')
nest.Connect(o, o, syn_spec={'synapse_model': 'rate_connection_delayed', 'receptor_type': 1})
t = nest.Create('step_rate_generator')
nest.Connect(t, o, syn_spec={'synapse_model': 'rate_connection_delayed', 'receptor_type': 2})
l = nest.Create('ctc_loss')
nest.Connect(o, [l.global_id] * len(o), conn_spec={'rule': 'one_to_one', 'make_symmetric': True},
             syn_spec={'synapse_model': 'gap_junction', 'receptor_type': range(1, len(o)+1)})
nest.Simulate(10)
