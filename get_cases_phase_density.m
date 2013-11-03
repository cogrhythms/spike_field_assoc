clear all; close all
%
% NOTE: the following paths are not defined. 
%       This is illustrative rather than working code.
%
addpath( '../../../cr' );
addpath( '../../../misc' );
base_nm   = 'data4cases';
n_cases   = 9;
dt        = 1e-3;   % seconds
%for c = 1 : n_cases
for c = [ 3 6 ]

  % =====================================================
  % Set the output name.
  % =====================================================
  out_nm        = sprintf( 'out/densities-case%d.mat', c );

  % =====================================================
  % For each case compute an estimate of the density
  % and save it.
  % =====================================================
  % Load the spiking.
  load_nm      = sprintf( 'mat/%s-case%d.mat', base_nm, c );
  load( load_nm     );

  % =====================================================
  % Load the lfp.
  % =====================================================
  load_nm      = sprintf( 'mat/%s-lfp.mat', base_nm );
  load( load_nm );

  % =====================================================
  % Get spike phases.
  % =====================================================
  f0        = 45; % Hz.
  bandwidth = 10; % Hz.
  b_lfp     = bandpass_lfp( lfp, dt, f0, bandwidth, ...
                            0, 'out/filter_chk.eps' );

  % =====================================================
  % Avoid transient of 100 samples at end of bandpassed
  % time-series.
  % =====================================================
  b_lfp     = b_lfp( 1:100, : );
  spikes_c1 = spikes_c1( 1:100, : );
  spikes_c2 = spikes_c2( 1:100, : );

  phases_c1 = get_spike_phases( spikes_c1, hilbert( b_lfp ));
  phases_c2 = get_spike_phases( spikes_c2, hilbert( b_lfp ));

  % =====================================================
  % Estimate the phase density for the 2 cases.
  % =====================================================
  n_domain_vals             = 400;    kappa = 20;
  [ density_est1 phase_domain1 ]  = cr_est_density( phases_c1(:), kappa, n_domain_vals  );
  [ density_est2 phase_domain2 ]  = cr_est_density( phases_c2(:), kappa, n_domain_vals  );

  nn = length( phases_c1 ); nn2 = length( phases_c2 );
  density_est1_ci = zeros( n_domain_vals/2, 30 );
  density_est2_ci = zeros( n_domain_vals/2, 30 );
  ii1             = unidrnd( nn, [ nn 30 ] );
  ii2             = unidrnd( nn2, [ nn2 30 ] );
  for jj = 1 : 30 
  
    ii = ii1( :, jj );
    [ density_est1_ci(:,jj) phase_domain1 ]  = cr_est_density( phases_c1(ii), kappa, n_domain_vals  );
    ii = ii2( :, jj );
    [ density_est2_ci(:,jj) phase_domain2 ]  = cr_est_density( phases_c2(ii), kappa, n_domain_vals  );

  end

  i_lo  = round( .025 * 30 );   i_hi = round( .975 * 30 );
  for n = 1 : n_domain_vals / 2

    svals                     = sort( density_est1_ci( n, : ), 'ascend' );
    density_est1_l95( n )     = svals( i_lo );
    density_est1_h95( n )     = svals( i_hi );

    svals                     = sort( density_est2_ci( n, : ), 'ascend' );
    density_est2_l95( n )     = svals( i_lo );
    density_est2_h95( n )     = svals( i_hi );
  end

  % =====================================================
  % save the results.
  % =====================================================
  fprintf( 'Saving to: %s\n', out_nm );
  save( out_nm, 'density_est1', 'phase_domain1', 'density_est2', 'phase_domain2', ...
                'density_est1_l95', 'density_est1_h95', ...
                'density_est2_l95', 'density_est2_h95' );


end
