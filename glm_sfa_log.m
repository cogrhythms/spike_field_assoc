% Publication:  A Procedure for Testing Across-Condition Rhythmic Spike-field Association Change
% Authors:      Kyle Q. Lepagea, Georgia G. Gregoriou, Mark A. Kramer, Mikio Aoi, Stephen J. Gotts, Uri T. Eden, Robert Desimone 
%
% Function:     glm_coh
% Input:
%
%     v     - cell array, each element another trial;
%             consisting of a nSamples x 2 dim. matrix.
%             first dimension is binned spike counts, 
%             second dimension is the lfp.
%     dt    - bin size for spikes and the sample period
%               for the lfp, in seconds.
%     bandwidth  - bandwidth to step by.
%
%     start_f - start frequency
%     stop_f  - stop frequency.
%     n_trim  - # of samples to remove from right hand side of data
%               sections to avoid filter edge effects.
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
%                         - dev ( 1 x no. of frequencies )
%                              The deviation, equal to -2 times the log likelihood.
%
%                         - stats ( 1 x no. of frequencies array of structures).  
%                             Each structure contains the stats structure
%                             returned by glm_fit_simple() at the corresponding frequency interval.
%
%
%                         - aic ( 1 x no. of frequencies )
%                             Aikiake Information Criterion = deviation  + 2 * no. of params.
%
%                         - f ( 1 x no. of frequencies )
%                             The centers of the frequency intervals of interest.  (Hz)
%
%                         - phases ( n_trials x n_bins_per_trial by no. of frequencies )  NOTE: When n_bins_per_trial is not 
%                                                                                               constant this dimension is a sum of
%                                                                                               the bins over trials.
%
%                             The instantaneous phase estimates as a function of frequency interval of interest.
%
%     
% 
function out = glm_sfa_log( v, dt, bandwidth, start_f, stop_f, n_trim, ...
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
  Ws          = 1.0/dt; 
  Wnyq        = Ws/2;           % Nyquist frequency, Hz

  % total_len unknown due to frequency dependent edge-effect trimming.
  filter_order  = 4;
  n_trim        = 100;  % Empirically determined in test_bp_filter_phase_edge_effect.m   
  if( abs( dt - 1e-3 ) > 1e-14 )
    fprintf( '\n\nCAUTION!!:  Filter edge effects tuned for bin size of 1 ms. n\n' )
    dbstack
    keyboard
  end


  spikes        = [];
  total_len     = sum( cell2mat( cellfun( @(x) size( x, 1 ), v, 'uniformoutput', false )));
  total_len     = total_len - n_trim * nTS;
  spikes        = zeros( total_len, 1 );
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
    if( 0 )
      Wlow                = startFreq; % Lower frequency bound = 4 Hz
      Whigh               = stopFreq; % Upper frequency bound is 8 Hz
      Wn                  = [Wlow/Wnyq+eps, Whigh/Wnyq-eps];  % Cutoff frequencies normalized to the
      [ bb, aa]            = butter(3,Wn);  % Design a 3rd order butterworth filter
      %[bb, aa]            = cheby2(10,20,Wn);  % Design a 3rd order butterworth filter
      [h{i_f},w]          = freqz(bb,aa);
      out.xfer_func{i_f}  = h{i_f};
      out.xfer_func_f     = w;
    else
 
      center_freq      = mean( [ startFreq stopFreq ] );
      Fs               = 1/dt;

      % NOTE: When specifying attenuation: filtfilt() is used which squares the magnitude response.
      %       This means that only half of the required attenuation needs to be specified.
      buffer           = 60;    % dB 
      fst1             = max( 0, startFreq - quarter_bw );       ast1  = ( 20*log10( range_lfp_dft ) + buffer ) / 2;   % dB
      fp1              = startFreq + quarter_bw;                 ap    = 2.5;                                          % dB
      fp2              = stopFreq  - quarter_bw;                 ast2  = ( 20*log10( range_lfp_dft ) + buffer ) / 2;   % dB
      fst2             = min( Wnyq, stopFreq + quarter_bw );
      filter_type      = 'ellip'
      try 
        if( fst1 < 1e-14 )

          fprintf( '\n\nThis code should never be executed.\n\n' );
          keyboard

          % IIR
          %spec_string      = 'n,fst1,fp1,fp2,fst2,ap';
          %ap_new           = 5;
          %bp_design        = fdesign.bandpass( spec_string, filter_order, fst1+eps, fp1, fp2, fst2, ap_new, Fs );
          %[ bb aa ]        = sos2tf( h_design.sosMatrix );

          % FIR
          spec_string      = 'N,Fc1,Fc2,Ast1,Ap,Ast2';
          bp_design        = fdesign.bandpass( spec_string, filter_order, fp1, fp2, ast1, ap, ast2, Fs );

          h_design         = design( bp_design, 'fir' );
          bb               = h_design.numerator;
          aa               = zeros( size( bb )); aa(1) = 1;

          [h{i_f},w]       = freqz(bb,aa);

        else

          if( strcmp( filter_type, 'iir' ))
              spec_string      = 'nb,na,fst1,fp1,fp2,fst2';
              bp_design        = fdesign.bandpass( spec_string, filter_order, filter_order, fst1, fp1, fp2, fst2, Fs );
              h_design         = design( bp_design, 'iir' );
              [ bb aa ]        = sos2tf( h_design.sosMatrix );
              [h{i_f},w]       = freqz(bb,aa);

          elseif( strcmp( filter_type, 'ellip' ))

              spec_string      = 'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2';
              bp_design        = fdesign.bandpass( spec_string, fst1, fp1, fp2, fst2, ast1, ap, ast2, Fs );
              h_design         = design( bp_design, 'ellip' );
              [ bb aa ]        = sos2tf( h_design.sosMatrix );
              [h{i_f},w]       = freqz(bb,aa);

          elseif( strcmp( filter_type, 'fir' ))
            % FIR
            spec_string      = 'N,Fc1,Fc2,Ast1,Ap,Ast2';
            bp_design        = fdesign.bandpass( spec_string, filter_order, fp1, fp2, ast1, ap, ast2, Fs );

            h_design         = design( bp_design, 'fir' );
            bb               = h_design.numerator;
            aa               = zeros( size( bb )); aa(1) = 1;

            [h{i_f},w]       = freqz(bb,aa);
          end
          out.xfer_func{i_f}  = h{i_f};
          out.xfer_func_f     = w;
        end
      catch me
        me.message
        dbstack
        keyboard
      end

      if( any( ~isfinite( h{i_f} )))
        fprintf( '\n\nPROBLEM: unstable bandpass filter.\n\n' );
        dbstack
        keyboard
      end
    end
    [ w_msg w_id ]  = lastwarn;
    if( ~isempty( w_msg ))
      keyboard
      [hh,w] = freqz(bb,aa);
      figure(1);clf;semilogy( w, abs(hh) ), print -depsc2 out/h.eps, close(1)
    end
    [hh,xfer_func_domain] = freqz(bb,aa);
    if( i_f == 1 )
      n_hh  = length( hh );
      bp_xfer_func = zeros( n_hh, nFreqs );
    end
    bp_xfer_func(:,i_f)   = hh;
    

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
      bp              = filtfilt( bb, aa, v{n}(:,2)-mean(v{n}(:,2)) );

      % =========================================================================
      % Trim the end of the bandpass filtered lfp to avoid the filter edge effect.
      % =========================================================================
      bp            = bp( 1 : N );


      % ==========================================
      % Some checks.
      % ==========================================
      [ w_msg w_id ]  = lastwarn;
      if( ~isempty( w_msg ))
        fprintf( sprintf( '\n\n%s\n\n', w_msg )) 
        keyboard
        [hh,w] = freqz(bb,aa);
        figure(1);clf;semilogy( w, abs(hh) ), print -depsc2 out/h.eps, close(1)
      end
      try 
        a            = hilbert( bp );
      catch me
        dbstack
        me.message
        keyboard
      end

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

    N           = total_len; cos_phase   = cos( phase ); sin_phase   = sin( phase );
    design_mat  = [ ones( N, 1 ) cos_phase sin_phase ];


    try
       [ b{i_f} dev{i_f} stats{i_f} ] = ...
                 glmfit( design_mat, spikes,   ...
                                'poisson', 'constant', 'off' );


    catch me
      dbstack
      me.message
      keyboard
    end
  end % for nFreqs

  if( debug_level == 2 )

    % =========================================
    % Testing for phase coherence across
    % the different frequencies.  If it
    % is present, across-frequency comparisons
    % will not be reliable.
    % =========================================
    nm_d      = sprintf( '%s.xfreq-phase-contamination.log-link.eps', f_out_in );
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
       %caxis( [ 0 2*pi*100*big_t(n_trials2use*nnn)/50 ] );colorbar
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

  end % debug_level2


  % ===================
  % Compile the output.
  % ===================
  out.coeffs              = cell2mat( b );  
  out.dev                 = cell2mat( dev );
  out.stats               = cell2mat( stats );
  out.aic                 = out.dev   +   2 * size( b{1}, 1 );
  out.f                   = f;  

  if( exist( 'phases' ))
    out.phases = phases;
  end

end % glm_sfa_log
