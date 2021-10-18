/* -*- c++ -*- */
/*                                  Apache License
 *                            Version 2.0, January 2004
 *                         http://www.apache.org/licenses/
 *
 *    TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION
 *
 *    1. Definitions.
 *
 *       "License" shall mean the terms and conditions for use, reproduction,
 *       and distribution as defined by Sections 1 through 9 of this document.
 *
 *       "Licensor" shall mean the copyright owner or entity authorized by
 *       the copyright owner that is granting the License.
 *
 *       "Legal Entity" shall mean the union of the acting entity and all
 *       other entities that control, are controlled by, or are under common
 *       control with that entity. For the purposes of this definition,
 *       "control" means (i) the power, direct or indirect, to cause the
 *       direction or management of such entity, whether by contract or
 *       otherwise, or (ii) ownership of fifty percent (50%) or more of the
 *       outstanding shares, or (iii) beneficial ownership of such entity.
 *
 *       "You" (or "Your") shall mean an individual or Legal Entity
 *       exercising permissions granted by this License.
 *
 *       "Source" form shall mean the preferred form for making modifications,
 *       including but not limited to software source code, documentation
 *       source, and configuration files.
 *
 *       "Object" form shall mean any form resulting from mechanical
 *       transformation or translation of a Source form, including but
 *       not limited to compiled object code, generated documentation,
 *       and conversions to other media types.
 *
 *       "Work" shall mean the work of authorship, whether in Source or
 *       Object form, made available under the License, as indicated by a
 *       copyright notice that is included in or attached to the work
 *       (an example is provided in the Appendix below).
 *
 *       "Derivative Works" shall mean any work, whether in Source or Object
 *       form, that is based on (or derived from) the Work and for which the
 *       editorial revisions, annotations, elaborations, or other modifications
 *       represent, as a whole, an original work of authorship. For the purposes
 *       of this License, Derivative Works shall not include works that remain
 *       separable from, or merely link (or bind by name) to the interfaces of,
 *       the Work and Derivative Works thereof.
 *
 *       "Contribution" shall mean any work of authorship, including
 *       the original version of the Work and any modifications or additions
 *       to that Work or Derivative Works thereof, that is intentionally
 *       submitted to Licensor for inclusion in the Work by the copyright owner
 *       or by an individual or Legal Entity authorized to submit on behalf of
 *       the copyright owner. For the purposes of this definition, "submitted"
 *       means any form of electronic, verbal, or written communication sent
 *       to the Licensor or its representatives, including but not limited to
 *       communication on electronic mailing lists, source code control systems,
 *       and issue tracking systems that are managed by, or on behalf of, the
 *       Licensor for the purpose of discussing and improving the Work, but
 *       excluding communication that is conspicuously marked or otherwise
 *       designated in writing by the copyright owner as "Not a Contribution."
 *
 *       "Contributor" shall mean Licensor and any individual or Legal Entity
 *       on behalf of whom a Contribution has been received by Licensor and
 *       subsequently incorporated within the Work.
 *
 *    2. Grant of Copyright License. Subject to the terms and conditions of
 *       this License, each Contributor hereby grants to You a perpetual,
 *       worldwide, non-exclusive, no-charge, royalty-free, irrevocable
 *       copyright license to reproduce, prepare Derivative Works of,
 *       publicly display, publicly perform, sublicense, and distribute the
 *       Work and such Derivative Works in Source or Object form.
 *
 *    3. Grant of Patent License. Subject to the terms and conditions of
 *       this License, each Contributor hereby grants to You a perpetual,
 *       worldwide, non-exclusive, no-charge, royalty-free, irrevocable
 *       (except as stated in this section) patent license to make, have made,
 *       use, offer to sell, sell, import, and otherwise transfer the Work,
 *       where such license applies only to those patent claims licensable
 *       by such Contributor that are necessarily infringed by their
 *       Contribution(s) alone or by combination of their Contribution(s)
 *       with the Work to which such Contribution(s) was submitted. If You
 *       institute patent litigation against any entity (including a
 *       cross-claim or counterclaim in a lawsuit) alleging that the Work
 *       or a Contribution incorporated within the Work constitutes direct
 *       or contributory patent infringement, then any patent licenses
 *       granted to You under this License for that Work shall terminate
 *       as of the date such litigation is filed.
 *
 *    4. Redistribution. You may reproduce and distribute copies of the
 *       Work or Derivative Works thereof in any medium, with or without
 *       modifications, and in Source or Object form, provided that You
 *       meet the following conditions:
 *
 *       (a) You must give any other recipients of the Work or
 *           Derivative Works a copy of this License; and
 *
 *       (b) You must cause any modified files to carry prominent notices
 *           stating that You changed the files; and
 *
 *       (c) You must retain, in the Source form of any Derivative Works
 *           that You distribute, all copyright, patent, trademark, and
 *           attribution notices from the Source form of the Work,
 *           excluding those notices that do not pertain to any part of
 *           the Derivative Works; and
 *
 *       (d) If the Work includes a "NOTICE" text file as part of its
 *           distribution, then any Derivative Works that You distribute must
 *           include a readable copy of the attribution notices contained
 *           within such NOTICE file, excluding those notices that do not
 *           pertain to any part of the Derivative Works, in at least one
 *           of the following places: within a NOTICE text file distributed
 *           as part of the Derivative Works; within the Source form or
 *           documentation, if provided along with the Derivative Works; or,
 *           within a display generated by the Derivative Works, if and
 *           wherever such third-party notices normally appear. The contents
 *           of the NOTICE file are for informational purposes only and
 *           do not modify the License. You may add Your own attribution
 *           notices within Derivative Works that You distribute, alongside
 *           or as an addendum to the NOTICE text from the Work, provided
 *           that such additional attribution notices cannot be construed
 *           as modifying the License.
 *
 *       You may add Your own copyright statement to Your modifications and
 *       may provide additional or different license terms and conditions
 *       for use, reproduction, or distribution of Your modifications, or
 *       for any such Derivative Works as a whole, provided Your use,
 *       reproduction, and distribution of the Work otherwise complies with
 *       the conditions stated in this License.
 *
 *    5. Submission of Contributions. Unless You explicitly state otherwise,
 *       any Contribution intentionally submitted for inclusion in the Work
 *       by You to the Licensor shall be under the terms and conditions of
 *       this License, without any additional terms or conditions.
 *       Notwithstanding the above, nothing herein shall supersede or modify
 *       the terms of any separate license agreement you may have executed
 *       with Licensor regarding such Contributions.
 *
 *    6. Trademarks. This License does not grant permission to use the trade
 *       names, trademarks, service marks, or product names of the Licensor,
 *       except as required for reasonable and customary use in describing the
 *       origin of the Work and reproducing the content of the NOTICE file.
 *
 *    7. Disclaimer of Warranty. Unless required by applicable law or
 *       agreed to in writing, Licensor provides the Work (and each
 *       Contributor provides its Contributions) on an "AS IS" BASIS,
 *       WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 *       implied, including, without limitation, any warranties or conditions
 *       of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
 *       PARTICULAR PURPOSE. You are solely responsible for determining the
 *       appropriateness of using or redistributing the Work and assume any
 *       risks associated with Your exercise of permissions under this License.
 *
 *    8. Limitation of Liability. In no event and under no legal theory,
 *       whether in tort (including negligence), contract, or otherwise,
 *       unless required by applicable law (such as deliberate and grossly
 *       negligent acts) or agreed to in writing, shall any Contributor be
 *       liable to You for damages, including any direct, indirect, special,
 *       incidental, or consequential damages of any character arising as a
 *       result of this License or out of the use or inability to use the
 *       Work (including but not limited to damages for loss of goodwill,
 *       work stoppage, computer failure or malfunction, or any and all
 *       other commercial damages or losses), even if such Contributor
 *       has been advised of the possibility of such damages.
 *
 *    9. Accepting Warranty or Additional Liability. While redistributing
 *       the Work or Derivative Works thereof, You may choose to offer,
 *       and charge a fee for, acceptance of support, warranty, indemnity,
 *       or other liability obligations and/or rights consistent with this
 *       License. However, in accepting such obligations, You may act only
 *       on Your own behalf and on Your sole responsibility, not on behalf
 *       of any other Contributor, and only if You agree to indemnify,
 *       defend, and hold each Contributor harmless for any liability
 *       incurred by, or claims asserted against, such Contributor by reason
 *       of your accepting any such warranty or additional liability.
 *
 *    END OF TERMS AND CONDITIONS
 *
 *    APPENDIX: How to apply the Apache License to your work.
 *
 *       To apply the Apache License to your work, attach the following
 *       boilerplate notice, with the fields enclosed by brackets "[]"
 *       replaced with your own identifying information. (Don't include
 *       the brackets!)  The text should be enclosed in the appropriate
 *       comment syntax for the file format. We also recommend that a
 *       file or class name and description of purpose be included on the
 *       same "printed page" as the copyright notice for easier
 *       identification within third-party archives.
 *
 *    Copyright [yyyy] [name of copyright owner]
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <gnuradio/math.h>
#include <volk/volk.h>
#include <boost/format.hpp>
#include <boost/math/special_functions/round.hpp>
#include <gnuradio/filter/pfb_arb_resampler.h>
#include <gnuradio/filter/firdes.h>
#include "corr_est_cc_impl.h"

namespace gr
{
  namespace ais
  {

    corr_est_cc::sptr
    corr_est_cc::make(const std::vector<gr_complex> &symbols,
                      float sps, unsigned int mark_delay,
                      float threshold)
    {
      return gnuradio::get_initial_sptr(new corr_est_cc_impl(symbols, sps, mark_delay, threshold));
    }

    /*
     * The private constructor
     */
    corr_est_cc_impl::corr_est_cc_impl(const std::vector<gr_complex> &symbols,
                                       float sps, unsigned int mark_delay,
                                       float threshold)
        : gr::sync_block("corr_est_cc",
                         gr::io_signature::make(1, 1, sizeof(gr_complex)),
                         gr::io_signature::make(1, 2, sizeof(gr_complex))),
          d_src_id(pmt::intern(alias()))
    {
      d_sps = sps;

      // Create time-reversed conjugate of symbols
      d_symbols = symbols;
      for (size_t i = 0; i < d_symbols.size(); i++)
      {
        d_symbols[i] = conj(d_symbols[i]);
      }
      std::reverse(d_symbols.begin(), d_symbols.end());

      d_mark_delay = mark_delay >= d_symbols.size() ? d_symbols.size() - 1
                                                    : mark_delay;

      // Compute a correlation threshold.
      // Compute the value of the discrete autocorrelation of the matched
      // filter with offset 0 (aka the autocorrelation peak).
      float corr = 0;
      for (size_t i = 0; i < d_symbols.size(); i++)
        corr += abs(d_symbols[i] * conj(d_symbols[i]));
      d_thresh = threshold * corr * corr;

      // Correlation filter
      d_filter = new kernel::fft_filter_ccc(1, d_symbols);

      // Per comments in gr-filter/include/gnuradio/filter/fft_filter.h,
      // set the block output multiple to the FFT filter kernel's internal,
      // assumed "nsamples", to ensure the scheduler always passes a
      // proper number of samples.
      int nsamples;
      nsamples = d_filter->set_taps(d_symbols);
      set_output_multiple(nsamples);

      // It looks like the kernel::fft_filter_ccc stashes a tail between
      // calls, so that contains our filtering history (I think).  The
      // fft_filter_ccc block (which calls the kernel::fft_filter_ccc) sets
      // the history to 1 (0 history items), so let's follow its lead.
      // set_history(1);

      // We'll (ab)use the history for our own purposes of tagging back in time.
      // Keep a history of the length of the sync word to delay for tagging.
      set_history(d_symbols.size() + 1);

      declare_sample_delay(1, 0);
      declare_sample_delay(0, d_symbols.size());

      // Setting the alignment multiple for volk causes problems with the
      // expected behavior of setting the output multiple for the FFT filter.
      // Don't set the alignment multiple.
      // const int alignment_multiple =
      //  volk_get_alignment() / sizeof(gr_complex);
      // set_alignment(std::max(1,alignment_multiple));

      // In order to easily support the optional second output,
      // don't deal with an unbounded max number of output items.
      // For the common case of not using the optional second output,
      // this ensures we optimally call the volk routines.
      const size_t nitems = 24 * 1024;
      set_max_noutput_items(nitems);
      d_corr = (gr_complex *)
          volk_malloc(sizeof(gr_complex) * nitems, volk_get_alignment());
      d_corr_mag = (float *)
          volk_malloc(sizeof(float) * nitems, volk_get_alignment());
    }

    /*
     * Our virtual destructor.
     */
    corr_est_cc_impl::~corr_est_cc_impl()
    {
      delete d_filter;
      volk_free(d_corr);
      volk_free(d_corr_mag);
    }

    std::vector<gr_complex>
    corr_est_cc_impl::symbols() const
    {
      return d_symbols;
    }

    void
    corr_est_cc_impl::set_symbols(const std::vector<gr_complex> &symbols)
    {
      gr::thread::scoped_lock lock(d_setlock);

      d_symbols = symbols;

      // Per comments in gr-filter/include/gnuradio/filter/fft_filter.h,
      // set the block output multiple to the FFT filter kernel's internal,
      // assumed "nsamples", to ensure the scheduler always passes a
      // proper number of samples.
      int nsamples;
      nsamples = d_filter->set_taps(d_symbols);
      set_output_multiple(nsamples);

      // It looks like the kernel::fft_filter_ccc stashes a tail between
      // calls, so that contains our filtering history (I think).  The
      // fft_filter_ccc block (which calls the kernel::fft_filter_ccc) sets
      // the history to 1 (0 history items), so let's follow its lead.
      // set_history(1);

      // We'll (ab)use the history for our own purposes of tagging back in time.
      // Keep a history of the length of the sync word to delay for tagging.
      set_history(d_symbols.size() + 1);

      declare_sample_delay(1, 0);
      declare_sample_delay(0, d_symbols.size());

      d_mark_delay = d_mark_delay >= d_symbols.size() ? d_symbols.size() - 1
                                                      : d_mark_delay;
    }

    int
    corr_est_cc_impl::work(int noutput_items,
                           gr_vector_const_void_star &input_items,
                           gr_vector_void_star &output_items)
    {
      gr::thread::scoped_lock lock(d_setlock);

      const gr_complex *in = (gr_complex *)input_items[0];
      gr_complex *out = (gr_complex *)output_items[0];
      gr_complex *corr;
      if (output_items.size() > 1)
        corr = (gr_complex *)output_items[1];
      else
        corr = d_corr;

      // Our correlation filter length
      unsigned int hist_len = history() - 1;

      // Delay the output by our correlation filter length so we can
      // tag backwards in time
      memcpy(out, &in[0], sizeof(gr_complex) * noutput_items);

      // Calculate the correlation of the non-delayed input with the
      // known symbols.
      d_filter->filter(noutput_items, &in[hist_len], corr);

      // Find the magnitude squared of the correlation
      volk_32fc_magnitude_squared_32f(&d_corr_mag[0], corr, noutput_items);

      int isps = (int)(d_sps + 0.5f);
      int i = 0;
      while (i < noutput_items)
      {
        // Look for the correlator output to cross the threshold
        if (d_corr_mag[i] <= d_thresh)
        {
          i++;
          continue;
        }
        // Go to (just past) the current correlator output peak
        while ((i < (noutput_items - 1)) &&
               (d_corr_mag[i] < d_corr_mag[i + 1]))
          i++;

        // Delaying the primary signal output by the matched filter
        // length using history(), means that the the peak output of
        // the matched filter aligns with the start of the desired
        // sync word in the primary signal output.  This corr_start
        // tag is not offset to another sample, so that downstream
        // data-aided blocks (like adaptive equalizers) know exactly
        // where the start of the correlated symbols are.
        add_item_tag(0, nitems_written(0) + i, pmt::intern("corr_start"),
                     pmt::from_double(d_corr_mag[i]), d_src_id);

        // Peak detector using a "center of mass" approach center
        // holds the +/- fraction of a sample index from the found
        // peak index to the estimated actual peak index.
        double center = 0.0;
        if (i > 0 and i < (noutput_items - 1))
        {
          double nom = 0, den = 0;
          for (int s = 0; s < 3; s++)
          {
            nom += (s + 1) * d_corr_mag[i + s - 1];
            den += d_corr_mag[i + s - 1];
          }
          center = nom / den - 2.0;
        }

        // Calculate the phase offset of the incoming signal.
        //
        // The analytic cross-correlation is:
        //
        // 2A*e_bb(t-t_d)*exp(-j*2*pi*f*(t-t_d) - j*phi_bb(t-t_d) - j*theta_c)
        //

        // The analytic auto-correlation's envelope, e_bb(), has its
        // peak at the "group delay" time, t = t_d.  The analytic
        // cross-correlation's center frequency phase shift, theta_c,
        // is determined from the argument of the analytic
        // cross-correlation at the "group delay" time, t = t_d.
        //
        // Taking the argument of the analytic cross-correlation at
        // any other time will include the baseband auto-correlation's
        // phase term, phi_bb(t-t_d), and a frequency dependent term
        // of the cross-correlation, which I don't believe maps simply
        // to expected symbol phase differences.
        float phase = fast_atan2f(corr[i].imag(), corr[i].real());
        int index = i + d_mark_delay;

        add_item_tag(0, nitems_written(0) + index, pmt::intern("phase_est"),
                     pmt::from_double(phase), d_src_id);
        add_item_tag(0, nitems_written(0) + index, pmt::intern("time_est"),
                     pmt::from_double(center), d_src_id);
        // N.B. the appropriate d_corr_mag[] index is "i", not "index".
        add_item_tag(0, nitems_written(0) + index, pmt::intern("corr_est"),
                     pmt::from_double(d_corr_mag[i]), d_src_id);

        if (output_items.size() > 1)
        {
          // N.B. these debug tags are not offset to avoid walking off out buf
          add_item_tag(1, nitems_written(0) + i, pmt::intern("phase_est"),
                       pmt::from_double(phase), d_src_id);
          add_item_tag(1, nitems_written(0) + i, pmt::intern("time_est"),
                       pmt::from_double(center), d_src_id);
          add_item_tag(1, nitems_written(0) + i, pmt::intern("corr_est"),
                       pmt::from_double(d_corr_mag[i]), d_src_id);
        }

        // Skip ahead to the next potential symbol peak
        // (for non-offset/interleaved symbols)
        i += isps;
      }

      // if (output_items.size() > 1)
      //   add_item_tag(1, nitems_written(0) + noutput_items - 1,
      //                pmt::intern("ce_eow"), pmt::from_uint64(noutput_items),
      //                d_src_id);

      return noutput_items;
    }

  } /* namespace ais */
} /* namespace gr */
