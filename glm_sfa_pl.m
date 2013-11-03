% Publication:  A Procedure for Testing Across-Condition Rhythmic Spike-field Association Change
% Authors:      Kyle Q. Lepagea, Georgia G. Gregoriou, Mark A. Kramer, Mikio Aoi, Stephen J. Gotts, Uri T. Eden, Robert Desimone 
%
% Function:     glm_sfa_pl
% Purpose:      GLM spike-field association with PL function.  
% Input:
%
%     v     - cell array, each element another trial;
%             consisting of a nSamples x 2 dim. matrix.
%             first dimension is binned spike counts, 
%             second dimension is the lfp.
%     dt    - bin size for spikes and the sample period
%               for the lfp, in seconds.
%
%    bandwidth  - size of frequency analysis interval. (Hz)
% 
%     start_f   - start of interval (Hz).
%     stop_f    - end of interval (Hz).
%     n_trim    - How many samples to discard at the end of every filtered data section to 
%                 avoid phase estimate bias due to filtering edge-effect.
% 
%     f_out_in      - name of file to write debugging and progress messages to.
%     debug_level   - determines the quantity of debug messages and plots written to file
%                     and/or displayed.  Set to zero to turn off debugging information.
%
% Output: 
%
%     out                 - structure with the following fields:
%     
%
%                         - xfer_func 
%                             Cell array consisting of an array for each analysis frequency interval.  Each array gives the filter
%                             response used to estimate the instantaneous phase.
%
%                         - xfer_func_f
%                             Cell array consisting of an array for each analysis frequency interval.  Each array gives the domain
%                             corresponding to the filter response in the corresponding field of out.xfer_func
%
%                         - coeffs  ( 3 x no. of frequencies )
%                             Matrix consisting of the beta_0, beta_s and beta_c coefficients, as a function of frequency.
%                             Frequency changes along the rows.
%
%                         - cov_mat (3 x 3 x no. of frequencies )
%                             3rd dimension specifies the analysis frequency interval.  The estimated covariance matrix of the
%                             beta coefficient estimators is given in the first 2 dimensions.  The coefficient ordering is 
%                             beta_0, beta_s and then beta_c so that the covariance between beta_s and beta_c estimators is 
%                             out.cov_mat( 3, 2, index_of_frequency_interval_of_interest )
%                             
%                         - sum_sqr_residuals ( 1 x no. of frequencies )
%                             Sum of the squared deviations: sum( ( counts - estimated_conditional_intensity / bin_width ).^2 ) as 
%                             a function of frequency interval of interest.
%
%                         - log_lik ( 1 x no of frequencies )
%                             Log of the likelihood.
%
%                         - n_params 
%                             no. of parameters.  Currently this number is 3, but in general other covariates could be added to
%                             the analysis.
%
%                         - aic ( 1 x no. of frequencies )
%                             Aikiake Information Criterion (AIC) = -2 * log( liklihood ) + 2 * n_params
%
%                         - score_eqn_chk_sum
%                             Sum of the absolute score equations.  This quantity is zero when the MLE is obtained.  In practice,
%                             this does not occur, but rather approximate MLE's are provided.  Lower values are preferred.
%
%                         - models ( n_trials x n_bins_per_trial by no. of frequencies )  NOTE: When n_bins_per_trial is not 
%                                                                                               constant this dimension is a sum of
%                                                                                               the bins over trials.
%
%                             The estimated conditional intensity of the n_param models to the spiking data as a function of 
%                             frequeny interval of interest.  
%
%                         - residuals ( n_trials x n_bins_per_trial by no. of frequencies )  NOTE: When n_bins_per_trial is not 
%                                                                                               constant this dimension is a sum of
%                                                                                               the bins over trials.
%
%                             As described in out.sum_sqr_residuals above. 
%
%                         - phases ( n_trials x n_bins_per_trial by no. of frequencies )  NOTE: When n_bins_per_trial is not 
%                                                                                               constant this dimension is a sum of
%                                                                                               the bins over trials.
%
%                             The instantaneous phase estimates as a function of frequency interval of interest.
%
%                         - f ( 1 x no. of frequencies )
%                             Center of the frequency intervals of interest.
%
%                         - modulation ( 1 x no. of frequencies )
%                             sqrt( beta_s^2 + beta_c^2 )
%
%                         - phase_rltn_w_field_deg ( 1 x no. of frequencies )
%                             atan2( beta_s, beta_c )
%
function out = glm_sfa_pl( v, dt, bandwidth, ...
                             start_f, stop_f, ...
                             n_trim, ...
                             f_out_in, debug_level )

  % Init.
  lastwarn( '' );
  if( nargin == 6 )
    f_out       = 1;
    debug_level = 0;
  elseif( nargin == 7 )
    f_out         = fopen( f_out_in, 'w+' );
    debug_level   = 0;
  elseif( nargin == 8 )
    f_out         = fopen( f_out_in, 'w+' );
  end

  nTS         = length( v );
  f           = start_f : bandwidth : stop_f - bandwidth;
  nFreqs      = length( f );
  Ws          = 1.0/dt; % Hz 
  Wnyq        = Ws/2;   % Nyquist frequency is Ws/2

  % total_len unknown due to frequency dependent edge-effect trimming.
  filter_order  = 4;
  %n_trim        = 100;  % Determine this empirically.
  fprintf( '\n =========================================================\n' );
  fprintf( 'n_trim set to %d.\n', n_trim );
  fprintf( '=========================================================\n\n' );
  if( abs( dt - 1e-3 ) > 1e-14 )
    fprintf( '\n\nCAUTION!!:  Filter edge effects tuned for bin size of 1 ms. n\n' )
    dbstack
    keyboard
  end

  % ==============================================
  % Find the dynamic range of the transformed LFP.
  % This is used to get the required attenuation
  % of the bandpass filter.
  % ==============================================
  total_len     = sum( cell2mat( cellfun( @(x) size( x, 1 ), v, 'uniformoutput', false )));
  total_len     = total_len - n_trim * nTS;
  spikes        = zeros( total_len, 1 );
  model         = zeros( total_len, nFreqs );
  residuals     = zeros( total_len, nFreqs );
  i_last_stop   = 0;
  max_dft       = -1;
  min_dft       = 1e10;
  for n = 1 : nTS
    N             = size( v{n}, 1 ) - n_trim;
    nnn           = N;

    i_start       = i_last_stop + 1;
    i_stop        = i_start + N - 1;
    i_last_stop   = i_stop;
    inds          = i_start:i_stop;
    spikes(inds)  = v{n}(1:N,1);

    dft_lfp       = fft( v{n}(:,2));
    max_dft       = max( max_dft, max( abs( dft_lfp )));
    min_dft       = min( min_dft, min( abs( dft_lfp )));
  end % for nTS
  range_lfp_dft   = max_dft - min_dft;
  fprintf( '\n\nRange LFP Mag. DFT: %.1f\n\n', range_lfp_dft );

  % ======================================
  % For each frequency band.
  % ======================================
  half_bandwidth  = bandwidth / 2;
  quarter_bw      = bandwidth / 4;
  for( i_f =  1 : nFreqs )

    % Estimate phase.
    phases{ i_f }       = [];
    startFreq           = f( i_f )  - half_bandwidth; % Hz
    stopFreq            = startFreq +      bandwidth; % Hz.
    try
    if( 0 ) % If a different filter is desired for estimating instantaneous phase.
      Wlow                = startFreq; % Lower frequency bound = 4 Hz
      Whigh               = stopFreq; % Upper frequency bound is 8 Hz
      Wn                  = [Wlow/Wnyq+eps, Whigh/Wnyq-eps];  % Cutoff frequencies normalized to the
      [bb, aa]            = butter(3,Wn);  % Design a 3rd order butterworth filter
      %[bb, aa]            = cheby2(2,20,Wn);  % Design a 3rd order butterworth filter
      [h{i_f},w]          = freqz(bb,aa);
      out.xfer_func{i_f}  = h{i_f};
      out.xfer_func_f     = w;
    else
      
      center_freq      = mean( [ startFreq stopFreq ] );
      Fs               = 1/dt;
      spec_string      = 'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2';
      %spec_string      = 'nb,na,fst1,fp1,fp2,fst2';

      % NOTE: When specifying attenuation: filtfilt() is used which squares the magnitude response.
      %       This means that only half of the required attenuation needs to be specified.
      buffer           = 60;    % dB 
      fst1             = max( 0, startFreq - quarter_bw );       ast1  = ( 20*log10( range_lfp_dft ) + buffer ) / 2;   % dB
      fp1              = startFreq + quarter_bw;                 ap    = 2.5;                                          % dB
      fp2              = stopFreq  - quarter_bw;                 ast2  = ( 20*log10( range_lfp_dft ) + buffer ) / 2;   % dB
      fst2             = min( Wnyq, stopFreq + quarter_bw );
       %bp_design        = fdesign.bandpass( spec_string, filter_order, fst1, fp1, fp2, fst2, ast1, ap, ast2, Fs );
       bp_design        = fdesign.bandpass( spec_string, fst1, fp1, fp2, fst2, ast1, ap, ast2, Fs );
       h_design         = design( bp_design, 'ellip' );
      %bp_design        = fdesign.bandpass( spec_string, filter_order, filter_order, fst1, fp1, fp2, fst2, Fs );
      %h_design         = design( bp_design, 'iir' );
      [ bb aa ]        = sos2tf( h_design.sosMatrix );
      [h{i_f},w]       = freqz(bb,aa);
      out.xfer_func{i_f}  = h{i_f};
      out.xfer_func_f     = w;
      if( any( ~isfinite( h{i_f} )))
        fprintf( '\n\nPROBLEM: unstable bandpass filter.\n\n' );
        dbstack
        keyboard
      end
    end
    catch me
      dbstack
      me
      keyboard
    end
    [ w_msg w_id ]  = lastwarn;
    if( ~isempty( w_msg ))
       keyboard
      [hhh,w] = freqz(bb,aa);
      figure(1);clf;semilogy( w, abs(hhh) ), print -depsc2 out/h.eps, close(1)
    end

    % Determine the total length of all of the time-series.
    phase               = zeros( total_len, 1 );
    i_last_stop         = 0;
    for( n = 1 : nTS )

      N             = length( v{n}( 1 : end - n_trim, 1 ));
      i_start       = i_last_stop + 1;
      i_stop        = i_start + N - 1;
      i_last_stop   = i_stop;
      inds          = i_start:i_stop;


      % ==========================================
      % Band-pass filter.
      % ==========================================
      bp            = filtfilt( bb, aa, v{n}(:,2)-mean(v{n}(:,2)) ); 

      % =========================================================================
      % Trim the end of the bandpass filtered lfp to avoid the filter edge effect.
      % =========================================================================
      bp            = bp( 1 : N );


      % ==========================================
      % Some checks.
      % ==========================================
      [ w_msg w_id ]  = lastwarn;
      if( ~isempty( w_msg ))
        fprintf( sprintf( '\n\nglm_sfa_pl(): %s\n\n', w_msg )) 
        [hhh,w] = freqz(bb,aa);
        figure(1);clf;semilogy( w, abs(hhh) ), print -depsc2 out/h.eps, close(1)
        keyboard
      end
      a     = hilbert( bp );
      if( any( isnan( bp )))  % Check if filter has gone awry.
        keyboard
      end

      phase(inds)  = angle( a );

    end % nTS


    if( debug_level >= 1 )
      % ==============================================
      % Phase estimated at adjacent frequencies can
      % be highly influenced by the lfp phase at other
      % frequencies due to "filter leakage".
      % Check this.
      % ==============================================
      phases{i_f} = phase;
    end


    N               = total_len;

    c_tapers_big    = cos( phase );
    s_tapers_big    = sin( phase );
    design_mat      = [ ones( N, 1 ) c_tapers_big s_tapers_big ];

    % ============================================================
    % Estimate the parameters in the model (ignoring the effect of 
    % the link function).
    % This is to prime the iterations with a reasonable first
    % guess at the parameter values.
    % ============================================================
    a               = sum( spikes ) / N;
    g_c               = c_tapers_big' * c_tapers_big;
    ig_c              = inv( g_c );
    g_s               = s_tapers_big' * s_tapers_big;
    ig_s              = inv( g_c );
    gamma_r           = ig_c *  c_tapers_big' * spikes;   % length K vector
    gamma_i           = ig_s *  s_tapers_big' * spikes;   % length K vector

    % Get the fit at the estimated quantities. 
    model(:,i_f)      = m_get_model( c_tapers_big, s_tapers_big, ...
                                      a, gamma_r, gamma_i, N, dt );

    % Get the score.
    score             = m_get_score( spikes, model(:,i_f), N, dt, c_tapers_big, s_tapers_big );
    score_array       = abs( [ score.bg score.c score.s ] );    % This will be near zero for MLE.
    try

      f0              = mean( [ startFreq stopFreq ] );
      fprintf( f_out, '\n\nglm_sfa_pl(): Refining parameter estimates for frequency %.1f.\n\n', f0 );
      fprintf( '\n\nglm_sfa_pl(): Refining parameter estimates for frequency %.1f.\n', f0 );

      % Refine the model parameter estimates.
      old_gamma_r   = gamma_r;  old_gamma_i = gamma_i;
      [ a gamma_r gamma_i lambdas hessian model(:,i_f) ] = m_refine_params( a, gamma_r, gamma_i, ...
                                                                            c_tapers_big, s_tapers_big, ...
                                                                            spikes, dt, N, model(:,i_f), f_out );


    catch me
      me.message
      keyboard
    end

    % =====================================================================
    % The score equations (3 of them, one for each parameter), are zero
    % at the MLE estimate and near zero near the MLE.
    % =====================================================================
    score                   = m_get_score( spikes, model(:,i_f), N, dt, c_tapers_big, s_tapers_big );


    % ================================================
    % Assign the fields in the output structure, out.
    % ================================================
    out.coeffs(1,i_f) = a;
    out.coeffs(2,i_f) = gamma_r;
    out.coeffs(3,i_f) = gamma_i;
    out.cov_mat(:,:,i_f) = -inv( hessian );
    residuals( :, i_f )  = spikes - dt * model(:,i_f);

    out.sum_sqr_residuals(i_f) = residuals(:,i_f)' * residuals(:,i_f);

    out.log_lik( i_f )            = -dt * sum( model(:,i_f) ) + spikes' * log( dt * model(:,i_f) );
    out.n_params                  = 3; %1 + 2 * length( gamma_r );
    out.aic( i_f )                = 2 * out.n_params - 2 * out.log_lik( i_f );
    out.score_eqn_chk_sum( i_f )  = abs( score.bg ) + abs( score.c ) + abs( score.s );

  end % for nFreqs
  out.models    = model;
  out.residuals = residuals;
  if( exist( 'phases' ))
    out.phases = phases;
  end

  if( debug_level == 2 )
    % =========================================
    % Testing for phase coherence across
    % the different frequencies.  If it
    % is present, across-frequency comparisons
    % will not be reliable.
    % =========================================
    nm_d      = sprintf( '%s.xfreq-phase-contamination.eps', f_out_in );
    phase_mat = cell2mat( phases );
    big_t     = [ 0:size(phase_mat,1) - 1 ]' *dt;
    unwrap_p  = unwrap( phase_mat );
    pred_unwrap_p = big_t * 2 * pi * f;
    n_trials2use  = 100 ; %min( [ 5 nTS ] );
    time          = n_trials2use * nnn * dt;
    figure(1);clf;
       subplot(1,3,1),imagesc( f, big_t, unwrap_p ),
       caxis( [ 0 2*pi*100*big_t(n_trials2use*nnn) ] );colorbar
       set(gca,'ylim',[0 time] )
       title( [{'Estimated LFP Phase vs'},{'LFP Rhythm Frequency'}] );
       ylabel( 'Time (in seconds, concatenated across trials )' )
       xlabel( 'Frequency (Hz)' )
       subplot(1,3,2),imagesc( f, big_t, unwrap_p - pred_unwrap_p ),
       caxis( [ -10  10 ] ); colorbar
       set(gca,'ylim',[0 time] )
       title( [{'Estimated LFP Phase -'}, {'Predicted LFP Phase'} ] );
       xlabel( 'Frequency (Hz)' )
       subplot(1,3,3),imagesc( f, big_t, diff( unwrap_p - pred_unwrap_p ) / ( 2 * pi * dt ) ),
       set(gca,'ylim',[0 time] )
       caxis( [ -20 20 ] ), colorbar
       title( 'Instantaneous Departure Frequency' );
       xlabel( 'Frequency (Hz)' )
       suptitle_cb_auto( 'BP Filter Phase Estimation Performance' )
       print( nm_d, '-depsc2' ); close(1);

  end

  % ===================
  % Compile the output.
  % ===================
  out.f                   = f;

  % ===================
  % Compute pvalues.
  % ===================
  for  jj = 1 : length( f )
    out.modulation( jj )            = sqrt( out.coeffs(2,jj)^2 + out.coeffs(3,jj)^2 ); 
    out.phase_rltn_w_field_deg(jj)  = atan2( out.coeffs(3,jj), out.coeffs(2,jj)) * 180/pi;
    if( debug_level >= 1 )
      fprintf( f_out, '\n\nglm_sfa_pl(): Frequency %.1f\t\tSum sqr. residuals = %.2e\n\n', ...
                          f( jj ), out.sum_sqr_residuals( jj ) ); 
      fprintf( f_out, '\n\nglm_sfa_pl(): n_params = %d, log_lik = %.2e, aic = %.2e\n\n', ...
                        out.n_params, out.log_lik( jj ), out.aic( jj ));
    end
  end

  if( nargin > 6 )
    fclose( f_out );
  end
 
  % Plots.
  if( debug_level == 2 )
    [ dummy i_max_freq ]    = max( out.modulation );
    n =  5;
    off = n*1000;
    figure(1);clf;
    plot( 50 * spikes, '-g', 'linewidth', 2 ); hold on
    plot( model(:,i_max_freq) ); hold on
    try
    plot( cell2mat( make_col( cellfun( @(x) x(:,2), v, 'uniformoutput', false )) ) + a/dt, '-r', ...
                    'linewidth', .5 ); hold on
    catch me
      keyboard
    end
    vline( [ 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 ], '-k' )
    set( gca, 'xlim', [ 1000 1200 ]+off )
    %set( gca, 'xlim', [  400 3600 ]+off )
    try
    print -depsc2 out/chk.eps, close(1)
    catch me 
      keyboard
    end
  end
end % glm_sfa_pl


function model  = m_get_model( c_tapers_big, s_tapers_big, a, gamma_r, gamma_i, N, dt )

    model = a + sum( c_tapers_big * gamma_r + s_tapers_big * gamma_i, 2 );
    model = model / dt;   % in Hz.
    model = max( model, eps );

end

% Refine parameter estimates using Newton-Raphson iteration.
% Problem is convex down (as required).  Thus, constrained extremum
% occurs on boundary, if constraint is not met.
function [ a gamma_r gamma_i lambdas hessian model ] = m_refine_params( a, gamma_r, gamma_i, ...
                                                                        c_tapers_big, s_tapers_big, ...
                                                                        spikes, dt, N, model, f_out )

  % ==========================================
  % log-lik undefined on constraint boundary.
  % Stay inside domain, and epsilon away
  % from boundary.
  % ==========================================
  % Just use eps, log( eps ) = -36.

  % ==========================================
  % If the model is anywhere negative, increase
  % the offset until it no longer is.
  % ==========================================
  bnd_proximity   = 1e-5;
  i_neg_spike     = find( model .*  spikes < 0.0 );
  i_neg_no_spike  = find( model .* ~spikes < 0.0 );
  min_model_spike = min( model( i_neg_spike ));
  if( ~isempty( i_neg_spike ))
    a               = a - min_model_spike*dt + eps;
    model           = m_get_model( c_tapers_big, s_tapers_big, ...
                                    a, gamma_r, gamma_i, N, dt );
    i_neg_spike           = find( model.*spikes < 0.0 );
  end
  lambdas         = [];
  tmp_lambdas     = [];
  assert( isempty( i_neg_spike ));

  i_spikes      = find( spikes > 0 );
  i_no_spikes   = find( spikes == 0 );
  i_too_small   = find( model( i_spikes ) < bnd_proximity );
  model( i_spikes( i_too_small ))   = bnd_proximity;


  % ===========================================
  % Perform unconstrained max. lik. estimation,
  % with Newton-Raphson using the Fourier sol'n.
  % as an initial guess.  Stay away from 
  % boundaries.
  % ===========================================
  n_max_iters = 7500;
  n_gamma = length( gamma_r );
  assert( n_gamma == length( gamma_i ) );
  chg               = zeros( 1 + 2 * n_gamma, n_max_iters );
  err               = zeros( 1 + 2 * n_gamma, n_max_iters );
  gc_scale          = 1;
  max_tapers        = max( c_tapers_big );
  n_params          = 1 + 2 * n_gamma;
  jj                = 0;
  b_continue        = true;
  tic
  while( b_continue )
    jj        = jj + 1;

    model     = m_get_model( c_tapers_big, s_tapers_big, ...
                             a, gamma_r, gamma_i, N, dt );

    model( i_spikes( i_too_small ))   = bnd_proximity;

    % =====================================
    % Don't use constraints.
    % =====================================
    grad_constraints  = [];
    c_err{jj}         = m_get_err( spikes, model, dt, N, c_tapers_big, s_tapers_big );

    c_mse( jj ) = c_err{jj}' * c_err{jj};
    hessian     = m_get_hessian( spikes, model, dt, N, ...
                                 c_tapers_big, s_tapers_big, ...
                                 a, gamma_r, gamma_i );

    % Hessian is added to derivatives wrt the gradient of the constraint.
    % But these are zero.
    n_constraints         = 0; %length( i_neg );
    J( 1:n_params, 1:n_params )                                                         = hessian;
    J( 1:n_params, n_params + 1 : n_params + n_constraints )                            = grad_constraints;
    J( n_params+1 : n_params + n_constraints, 1:n_params )                              = grad_constraints';
    J( n_params+1 : n_params + n_constraints, n_params+1 : n_params + n_constraints )   = zeros( n_constraints );

    [ ju jv ] = eig( J );
    dd        = diag( jv );
    if( max( dd ) / min( dd ) > 1e15 )
      fprintf( f_out, 'CAUTION: Newton Raphson about to take a big step.\n' );
      keyboard
    end

    iJ                            = inv( J );
    c_step{jj}                    = -iJ * c_err{jj};
    c_chg{jj}                     = [ a; gamma_r; gamma_i; lambdas ] + c_step{jj};
    tmp_a                         = c_chg{jj}( 1 );    
    tmp_gamma_r                   = c_chg{jj}( 2:n_gamma+1 );
    tmp_gamma_i                   = c_chg{jj}( n_gamma+2:n_params );
    lambdas                       = c_chg{jj}( n_params+1:end );

    model         = m_get_model( c_tapers_big, s_tapers_big, ...
                                  tmp_a, tmp_gamma_r, tmp_gamma_i, N, dt );

    model( i_spikes( i_too_small ))   = bnd_proximity;

    % Find where model hits boundaries.
    % Back-off until all constraints valid.
    % Remember where the constraints were violated, and enforce the
    % gradient of the log-lik and the linear combo of constraints
    % to be parallel. 
    %i_new_neg       = make_col( setdiff( find( model .* spikes + eps < eps ), i_neg ));
    %i_neg           = [ i_neg; i_new_neg ];
    while( min( model .* spikes ) < 0.0 )
      fprintf( '\n\nShould not happen!!!\n\n' );
      keyboard
      c_step{jj}    = c_step{jj} / 2;
      c_chg{jj}     = [ a; gamma_r; gamma_i ; lambdas ] + c_step{jj};
      tmp_a         = c_chg{jj}( 1 );
      tmp_gamma_r   = c_chg{jj}( 2 : n_gamma + 1 );
      tmp_gamma_i   = c_chg{jj}( n_gamma + 2 : n_params );
      tmp_lambdas   = c_chg{jj}( n_params+1 : end );
      model         = m_get_model( c_tapers_big, s_tapers_big, ...
                                    tmp_a, tmp_gamma_r, tmp_gamma_i, N, dt );
    end
    a           = tmp_a;
    gamma_r     = tmp_gamma_r;
    gamma_i     = tmp_gamma_i;
    lambdas     = tmp_lambdas;

    gamma_rs(:,jj)  = gamma_r;
    gamma_is(:,jj)  = gamma_i;

    aa(jj)      = a;
    if( jj ~= 1 )
      if( c_mse(jj) < 1e-5 )
        b_continue = false;
      end
    end
    if( jj > n_max_iters - 1 )
      b_continue = false;
    end

    if( b_continue == false )
      fprintf( '\n\nglm_sfa_pl(): Final jj: %d\n', jj );
      fprintf( 'glm_sfa_pl(): Final sum squared error: %.3e\n', c_mse(jj) );
      if( jj == n_max_iters )
        fprintf( 'glm_sfa_pl(): maximum iterations attained.\n\n' );
      end
    end
  end % while  b_continue
  toc

  i_neg                         = find( model( i_no_spikes ) < 0 );
  model( i_no_spikes( i_neg ))  = eps;
  hessian   = m_get_hessian( spikes, model, dt, N, ...
                             c_tapers_big, s_tapers_big, ...
                             a, gamma_r, gamma_i );


  fprintf( f_out, '\n\nglm_sfa_pl(): Newton-Raphson MSE: before iterations %.2e, after iterations %.2e\n\n', c_mse(1), c_mse(end) );
end


function [ err grad_constraints ] = m_get_err_constrained( spikes, model, lambdas, ...
                                              dt, N, c_tapers_big, s_tapers_big, ...
                                              a, gamma_r, gamma_i, max_tapers, i_neg, ...
                                              gc_scale )

  i_pos   = find( model > eps ); n_pos = length( i_pos );
  err     = -dt * n_pos                                   + spikes'   * ( 1.0 ./ model );
  err     = [ err ( -dt * sum( c_tapers_big( i_pos, : ) ) + ( spikes(i_pos) ./ model(i_pos) )' * c_tapers_big(i_pos,:) ) ];
  err     = [ err ( -dt * sum( s_tapers_big( i_pos, : ) ) + ( spikes(i_pos) ./ model(i_pos) )' * s_tapers_big(i_pos,:) ) ]';

  n_params                        = 1 + 2 * length( gamma_r );

  % Never add constraints.  Solution lies within the domain. 
  grad_constraints                = [];

end


function hessian   = m_get_hessian( spikes, model, dt, N, ...
                                    c_tapers_big, s_tapers_big, ...
                                    a, gamma_r, gamma_i )
  % Init.
  n_gamma         = length( gamma_r );
  model           = model * dt;         % counts / bin.
  model_sqr       = model.^2;

  % d2l_da2
  d2l_da2         = -spikes' * model_sqr;

  % d_m_d_gamma_x_j :   derivative of model for each time index j as a function 
  %                     of parameter gamma_r and gamma_i.
  d_m_d_gamma_r_j     = c_tapers_big;
  d_m_d_gamma_i_j     = s_tapers_big;

  % New way.  Modified to incorporate characteristic function
  % on May 9, 2011.
  char_func           = model > 0;
  max_model_sqr_w_0   = max( model_sqr, 0 );
  sqrt_D      = sqrt( spikes ./ max_model_sqr_w_0 .* char_func );
  sqrt_D_mat  = sqrt_D * ones( 1, n_gamma );
  H           = [ ones(N,1) c_tapers_big s_tapers_big ];
  sqrt_D_mat2 = sqrt_D * ones( 1, 1 + 2 * n_gamma );
  hessian     = -( sqrt_D_mat2 .* H )' * ( sqrt_D_mat2 .* H );
  if( any( hessian(:)) > 1e30 )
    dbstack
    keyboard
  end
  %err_cov   = inv( -hessian );
end


function err = m_get_err( spikes, model, dt, N, c_tapers_big, s_tapers_big )
  i_nneg  = find( model > eps );
  err     = -dt * length( i_nneg )                      +   spikes( i_nneg )'                      * ( 1.0 ./ model( i_nneg ));
  err     = [ err ( -dt * sum( c_tapers_big( i_nneg ) ) + ( spikes( i_nneg ) ./ model( i_nneg ) )' * c_tapers_big( i_nneg ) ) ];
  err     = [ err ( -dt * sum( s_tapers_big( i_nneg ) ) + ( spikes( i_nneg ) ./ model( i_nneg ) )' * s_tapers_big( i_nneg ) ) ]';
end

function    score   = m_get_score( spikes, model, N, dt, c_tapers_big, s_tapers_big )
    i_nneg    = find( model > eps );
    score.bg  =   -dt * length(i_nneg)               +   spikes( i_nneg )' * ( 1.0 ./ model( i_nneg ) );
    score.c   = ( -dt * sum( c_tapers_big( i_nneg )) + ( spikes( i_nneg ) ./ model( i_nneg ) )' * c_tapers_big( i_nneg ) );
    score.s   = ( -dt * sum( s_tapers_big( i_nneg )) + ( spikes( i_nneg ) ./ model( i_nneg ) )' * s_tapers_big( i_nneg ) );
end
