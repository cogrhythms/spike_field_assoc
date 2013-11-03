% Publication:  A Procedure for Testing Across-Condition Rhythmic Spike-field Association Change
% Authors:      Kyle Q. Lepagea, Georgia G. Gregoriou, Mark A. Kramer, Mikio Aoi, Stephen J. Gotts, Uri T. Eden, Robert Desimone 
%
% Script:       get_cases_param_est_with_exp_many_f()
% Purpose:      An example analysis of the methodology applied to synthetic data. 
%
clear all; close all
if ispc
  system( 'del out\msc.ps' );
else
  system( 'rm out/msc.ps' );
end
base_nm   = 'data4cases';
n_cases   = 11;
dt        = 1e-3;   % seconds
for c = [ 3 6 ] %^1 : n_cases

  fprintf( '\n\nAnalyzing Case %d\n\n', c );

  % =====================================================
  % Set the output name.
  % =====================================================
  out_nm        = sprintf( 'out/est_params-case%d_many_f.mat', c );

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
  % Place the data into the v form.
  %
  % spikes_cx   is no. of bins by no. of trials.
  % lfp         is no. of samples by no. of trials.
  %
  % Note that the length of each trial must be the same
  % here, but that the v format allows for unequal length
  % trials. glm_sfa_pl() and glm_sfa_log() work with
  % unequal length trials.
  % =====================================================
  v1  = m_form_v( spikes_c1, lfp );
  v2  = m_form_v( spikes_c2, lfp );

  % =================================================================
  % Estimate the coupling at f0 with the piece-wise linear function.
  % =================================================================
  bandwidth   = 10; % Hz.
  start_f     = 10;     % Hz
  stop_f      = 100.0;    % Hz
  n_trim      = 100;
  f_out       = 'out/glm_hs_simple.out'; 
  debug_level = 1;  % To get phase output set >= 1.

   pl_out1     = glm_sfa_pl( v1, dt, ...
                                   bandwidth, start_f, stop_f, n_trim, ...
                                   f_out, debug_level );

   pl_out2     = glm_sfa_pl( v2, dt, ...
                                   bandwidth, start_f, stop_f, n_trim, ...
                                   f_out, debug_level );

   keyboard
  coeffs_1    = pl_out1.coeffs;      coeffs_2    = pl_out2.coeffs;
  cov_mat_1   = pl_out1.cov_mat;     cov_mat_2   = pl_out2.cov_mat;

  % =========================================================
  % Estimated the coupling with the log-link function.
  % =========================================================
  bandwidth   = 10; % Hz.
  start_f     = 10;     % Hz
  stop_f      = 100;    % Hz
  f_out       = 'out/glm_exp.out'; % Currently not written to for any
                                   % debug_level.  Used as a base
                                   % file name for debug_level=2 output.
  debug_level = 1;
  exp_out_1   = glm_sfa_log( v1, dt, bandwidth, ...
                                         start_f, stop_f, n_trim, ...
                                         f_out, debug_level );

  exp_out_2   = glm_sfa_log( v2, dt, bandwidth, ...
                                         start_f, stop_f, n_trim, ...
                                         f_out, debug_level );


    % =====================================================
    % Test for across condition modulation change.
    % Piece-wise linear model.
    % =====================================================
    n_domain_vals   = 500;
    tic
    for i_f = 1 : size( coeffs_1, 2 )

      pl_mod1( i_f )     = sqrt( coeffs_1(2:end, i_f)' * coeffs_1(2:end, i_f) );
      pl_mod2( i_f )     = sqrt( coeffs_2(2:end, i_f)' * coeffs_2(2:end, i_f) );
      pl_sigma1( i_f )   = sqrt( mean( diag( cov_mat_1(2:end,2:end,i_f) )) );
      pl_sigma2( i_f )   = sqrt( mean( diag( cov_mat_2(2:end,2:end,i_f) )) );
      phi1            = atan2( coeffs_1(3,i_f), coeffs_1(2,i_f) );
      phi2            = atan2( coeffs_2(3,i_f), coeffs_2(2,i_f) );
      [ p_val12_pl( i_f ) null_pdf null_pdf_x test_stat_pl( i_f ) ] = ...
          test_modulation_diff( pl_mod1( i_f ), pl_mod2( i_f ), pl_sigma1( i_f ), pl_sigma2( i_f ), n_domain_vals );

      % Check for normalization:
      norm_chk_pl_mod_test( i_f )        = sum( null_pdf ) * ( null_pdf_x(2) - null_pdf_x(1) );
      fprintf( 'Integral of null PDF: %.3e\n', norm_chk_pl_mod_test( i_f ));
    end
    toc

    % =====================================================
    % Test for across condition modulation change.
    % log linear model.
    % =====================================================
    n_domain_vals   = 500;
    for i_f = 1 : length( exp_out_1.f )

      e_mod1( i_f )     = sqrt( exp_out_1.coeffs(2:end, i_f)' * exp_out_1.coeffs(2:end, i_f) );
      e_mod2( i_f )     = sqrt( exp_out_2.coeffs(2:end, i_f)' * exp_out_2.coeffs(2:end, i_f) );
      cov_mat_1         = exp_out_1.stats( i_f ).covb;    cov_mat_2         = exp_out_2.stats( i_f ).covb;
      e_sigma1( i_f )   = sqrt( mean( diag( cov_mat_1(2:end,2:end) )) );
      e_sigma2( i_f )   = sqrt( mean( diag( cov_mat_2(2:end,2:end) )) );
      [ p_val12_log( i_f ) null_pdf null_pdf_x test_stat_log( i_f ) ] = ...
          test_modulation_diff( e_mod1( i_f ), e_mod2( i_f ), e_sigma1( i_f ), e_sigma2( i_f ), n_domain_vals );

      % Check for normalization:
      norm_chk_exp_mod_test( i_f )        = sum( null_pdf ) * ( null_pdf_x(2) - null_pdf_x(1) );
      %fprintf( '\n\nIntegral of null PDF: %.3e\n\n', norm_chk_exp_mod_test( i_f ));
    end


 
  fprintf( '\n\nSaving to %s.\n\n', out_nm );
  save( out_nm, 'exp_out_1', 'exp_out_2', 'pl_out1', 'pl_out2', ...
                  'p_val12_pl', 'test_stat_pl', 'p_val12_log', 'test_stat_log', ...
                  'norm_chk_pl_mod_test', 'norm_chk_exp_mod_test' );
  fprintf( '\n\nDone.\n\n' );
end


