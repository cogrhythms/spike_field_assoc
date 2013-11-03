function v = m_form_v( spikes, rate )
  n_trials = size( spikes, 2 );
  for n = 1 : n_trials
    v{n} = [ spikes(:,n) rate(:,n) ];
  end
end
