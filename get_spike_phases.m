% b_lfp   - analytic signal of bandpass filtered lfp.
% spikes  - spikes.
% Note: multiple spikes in a single bin results in a single
%       phase value.
function phases = get_spike_phases( spikes, ab_lfp )

  [ n n_trials ]    = size( spikes );
  [ nc n_trialsc ]  = size( ab_lfp );

  assert( n == nc && n_trials == n_trialsc );

  phases            = angle( ab_lfp( spikes > 0 ));
end
