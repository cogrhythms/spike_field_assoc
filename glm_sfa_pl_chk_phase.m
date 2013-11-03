% Function:     glm_sfa_pl_chk_phase 
% Purpose:      To explore the phase estimates vs. frequency.  
%               Quality control.  This is a companion function
%               to glm_sfa_pl() that provides further
%               debuging without the computational overhead.
% Input:
%
%     v       - Data used in call to glm_sfa_pl().
%
%     dt      - bin size for spikes and the sample period
%               for the lfp, in seconds.
%
%     out     - output structure returned from glm_sfa_pl().
%
%     n_secs  - duration of time to investigate the unwrapped phase.  
%               Note that unwrapping occurs across trials.
% 
% 
%     fname   - name of file to write debugging and progress messages to.
%
% Output: 
%
%     
%
function debug_out = glm_sfa_pl_chk_phase( v, dt, out, n_secs, fname )


  % Init.
  debug_out         = 1
  n_seconds         = n_secs;
  debug_level       = 2;
  lastwarn( '' );

  out_nm            = sprintf( '%s.phase-check.ps', fname );
  fprintf( '\n\n==========================================\n' );
  fprintf( 'Preparing to delete %s\n', out_nm );
  if ispc
    system( sprintf( 'del %s', out_nm ));
  else
    system( sprintf( 'rm %s', out_nm ));
  end
  fprintf( '\t %s deleted.\n', out_nm );
  fprintf( '==========================================\n\n' );

  nTS         = length( v );
  nFreqs      = length( out.f );
  Ws          = 1.0/dt; % Sampling frequency is 600.615 Hz
  Wnyq        = Ws/2; % Nyquist frequency is Ws/2
  total_len   = sum( cell2mat( cellfun( @(x) size( x, 1 ), v, 'uniformoutput', false )));
  spikes      = zeros( total_len, 1 );
  i_last_stop = 0;
  max_dft     = -1;
  min_dft     = 1e10;
  for n = 1 : nTS
    N             = size( v{n}, 1 );
    i_start       = i_last_stop + 1;
    i_stop        = i_start + N - 1;
    i_last_stop   = i_stop;
    inds          = i_start:i_stop;
    spikes(inds)  = v{n}(:,1);
    dft_lfp       = fft( v{n}(:,2) - mean( v{n}(:,2)));
    max_dft       = max( max_dft, max( abs( dft_lfp )));
    min_dft       = min( min_dft, min( abs( dft_lfp )));
  end % for nTS

  range_lfp_dft   = max_dft - min_dft;
  fprintf( '\n\nRange LFP Mag. DFT: %.1f\n\n', range_lfp_dft );

  % =========================================
  % Phase filter vs frequency. Mag. response.
  % =========================================
  h_img     = cell2mat( out.xfer_func );
  figure(1);clf;
    imagesc( out.f, out.xfer_func_f/pi*.5/dt, 20*log10( abs( h_img ))); colorbar;
    title( [{'Mag. Response (dB) vs. Frequ.'}, ...
            {sprintf('Worst-case Reqd. Buffer = %.1f dB', 20*log10( range_lfp_dft ))} ]  );
    xlabel( 'Frequency (Hz)' );
    ylabel( 'Frequency (Hz)' );
    set(gca,'ylim',[0 max(out.f)] );
    print( out_nm, '-dpsc2', '-append' ); close(1);

  % =========================================
  % =========================================
  phase_mat     = cell2mat( out.phases  );
  big_t         = [ 0:size(phase_mat,1) - 1 ]' *dt;
  unwrap_p      = unwrap( phase_mat );
  pred_unwrap_p = big_t * 2 * pi * out.f;

  figure(1);clf;
    subplot(1,3,1),imagesc( out.f, big_t, unwrap_p );
      caxis( [ 0 2*pi*100*big_t(round(n_seconds/dt)) ] );colorbar
      set(gca,'ylim',[0 n_seconds] );
      title( [{'Estimated LFP Phase vs'},{'LFP Rhythm Frequency'}] );
      ylabel( 'Time (in seconds, concatenated across trials )' );
      xlabel( 'Frequency (Hz)' );
    subplot(1,3,2),imagesc( out.f, big_t, unwrap_p - pred_unwrap_p );
      caxis( [ 0 2*pi*100*big_t(round(n_seconds/dt))/10 ] );colorbar;
      set(gca,'ylim',[0 n_seconds] );
      title( [{'Estimated LFP Phase -'}, {'Predicted LFP Phase'} ] );
      xlabel( 'Frequency (Hz)' );
    subplot(1,3,3),imagesc( out.f, big_t, diff( unwrap_p - pred_unwrap_p ) / ( 2 * pi * dt ) ),
      set(gca,'ylim',[0 n_seconds] );
      caxis( [ -20 20 ] ), colorbar;
      title( 'Instantaneous Departure Frequency' );
      xlabel( 'Frequency (Hz)' );
      suptitle_cb_auto( 'BP Filter Phase Estimation Performance' );
    print( out_nm, '-dpsc2', '-append' ); close(1);

  fprintf ('\n\n==================================\n' );
  fprintf ('Done.\n' );
  fprintf ('==================================\n\n' );

end
