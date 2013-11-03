% Test for difference in modulation.
% 
% a is sqrt( x^2 + y^2 ), X,Y normal standard dev. sigma_a.
% b is sqrt( u^2 + v^2 ), U,V normal standard dev. sigma_b.
% a and b are independent.
% 
% domain goes with pdf_a and pdf_b .
% domain2 goes with pdf_s
%
% s is the measured value (a - b) to be tested.
% it is the difference in modulations.
%
% When the variances differ dramatically, the 
% domain can become unreasonable for one of
% the two conditions.  Thus, to have the same
% discretized grid can involve a very large number
% of values, and at the extreme of the domain numerical
% problems can become apparent.  To avoid this, 
% use two different domains for the two different pdfs
% and then interpolate to the coarser for the computation
% of the null density.
%
function [ p_val pdf_s domain2 s ] = ...
            test_modulation_diff( a, b, sigma_a, sigma_b, n_domain_vals )


  % ==================================================
  % Get the distribution of  a - b, under the null
  % (under null, a = b).  
  % ==================================================
  m_sigma   = max( [ sigma_a sigma_b ] );
  upper     = max( [ a b ] ) + 10 * m_sigma;
  lower     = 0; 

  % ==================================================
  % statistic, s = a - b = 0 under null.
  % ==================================================
  s         = a - b;

  % ==================================================
  % Get the densities.
  % ==================================================
  domain    = linspace( lower, upper, n_domain_vals );
  pdf_a     = m_get_pdf( a, sigma_a, domain );
  pdf_b     = m_get_pdf( a, sigma_b, domain );

  % ==================================================================
  % Perform the non-flipped convolution (correlation-typ-computation).
  % Yields pdf of s given a = b.
  % ==================================================================
  domain2   = [ -fliplr( domain(2:end)) domain ];
  dx        = domain2(2)-domain2(1);

  % =========================================
  % Augment pdf_a with zeros for negative
  % domain.
  % =========================================
  n_z         = n_domain_vals - 1;
  pdf_a_full  = [ zeros(1,n_z) pdf_a ];
  conv_mat    = hankel( pdf_a_full );
  conv_mat    = conv_mat( :, 1:n_domain_vals );
  pdf_s       = dx * conv_mat * pdf_b';
  norm_chk  = sum( pdf_s ) * dx;


  % ==================================================================
  % Compute the test p-values.
  % ==================================================================
  p_val_lo  = sum( pdf_s( domain2 < -abs( s ))) * dx;
  p_val_hi  = sum( pdf_s( domain2 > abs( s ))) * dx;
  p_val     = p_val_lo + p_val_hi;

  if( 0 )
    figure(1);clf;
      subplot(1,1,1),
      plot( domain, pdf_a, '-b', domain, pdf_b, '-g', domain2, pdf_s, '-k' );
        vline( s, '-r' );
    %subplot(1,2,2),
      %imagesc( domain, domain2, conv_mat ), colorbar
      print -depsc2 out/pdf.eps, close(1)
  end % if( 0 )

  if( ~isfinite( norm_chk ))
    fprintf( '\n\ntest_modulation_diff(): Numerical issue detected.\n' );
    fprintf( 'Reporting maximum cantelli intersection bound for test p-value.\n\n' );

    %c_domain                = [ a b ]; %linspace( a, b, 500 );
    distance_from_a         = abs( b - a );
    distance_from_b         = distance_from_a;
    distance_from_a_in_std  = distance_from_a / sigma_a;
    distance_from_b_in_std  = distance_from_b / sigma_b;
    p_a_gt_c_upper_bnd      = 1.0 ./ ( 1.0 + distance_from_a_in_std.^2 );
    p_b_gt_c_upper_bnd      = 1.0 ./ ( 1.0 + distance_from_b_in_std.^2 );
    p_val                   = max( [ p_a_gt_c_upper_bnd p_b_gt_c_upper_bnd ] );
      
  end % terrible numerics
end % function

function pdf = m_get_pdf( v, sigma, domain )
    x     = domain;
    var_  = sigma^2;
    try
      pdf   = x/var_  .* exp( -( x.^2 + v^2 )/(2*var_ )) .* besseli(0,x*v/var_ );
    catch me
      dbstack
      me.message
      keyboard
    end
end
