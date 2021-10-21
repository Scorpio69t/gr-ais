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
#include "msk_timing_recovery_cc_impl.h"
#include <gnuradio/math.h>
#include <gnuradio/filter/firdes.h>

namespace gr
{
  namespace ais
  {

    msk_timing_recovery_cc::sptr
    msk_timing_recovery_cc::make(float sps, float gain, float limit, int osps = 1)
    {
      return gnuradio::get_initial_sptr(new msk_timing_recovery_cc_impl(sps, gain, limit, osps));
    }

    /*
     * The private constructor
     */
    msk_timing_recovery_cc_impl::msk_timing_recovery_cc_impl(float sps, float gain, float limit, int osps)
        : gr::block("msk_timing_recovery_cc",
                    gr::io_signature::make(1, 1, sizeof(gr_complex)),
                    gr::io_signature::make3(1, 3, sizeof(gr_complex), sizeof(float), sizeof(float))),
          d_limit(limit),
          d_interp(new filter::mmse_fir_interpolator_cc()),
          d_dly_conj_1(0),
          d_dly_conj_2(0),
          d_dly_diff_1(0),
          d_mu(0.5),
          d_div(0),
          d_osps(osps)
    {
      set_sps(sps);
      enable_update_rate(true); // fixes tag propagation through variable rate blox
      set_gain(gain);
      if (d_osps != 1 && d_osps != 2)
        throw std::out_of_range("osps must be 1 or 2");
    }

    /*
     * Our virtual destructor.
     */
    msk_timing_recovery_cc_impl::~msk_timing_recovery_cc_impl()
    {
      delete d_interp;
    }

    void msk_timing_recovery_cc_impl::set_sps(float sps)
    {
      d_sps = sps / 2.0; // loop runs at 2x sps
      d_omega = d_sps;
      set_relative_rate(d_osps / sps);
      //        set_history(d_sps);
    }

    float msk_timing_recovery_cc_impl::get_sps(void)
    {
      return d_sps;
    }

    void msk_timing_recovery_cc_impl::set_gain(float gain)
    {
      d_gain = gain;
      if (d_gain <= 0)
        throw std::out_of_range("Gain must be positive");
      d_gain_omega = d_gain * d_gain * 0.25;
    }

    float msk_timing_recovery_cc_impl::get_gain(void)
    {
      return d_gain;
    }

    void msk_timing_recovery_cc_impl::set_limit(float limit)
    {
      d_limit = limit;
    }

    float msk_timing_recovery_cc_impl::get_limit(void)
    {
      return d_limit;
    }

    void
    msk_timing_recovery_cc_impl::forecast(int noutput_items, gr_vector_int &ninput_items_required)
    {
      /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
      unsigned ninputs = ninput_items_required.size();
      for (unsigned i = 0; i < ninputs; i++)
      {
        ninput_items_required[i] = (int)ceil((noutput_items * d_sps * 2) + 3.0 * d_sps + d_interp->ntaps());
      }
    }

    int
    msk_timing_recovery_cc_impl::general_work(int noutput_items,
                                              gr_vector_int &ninput_items,
                                              gr_vector_const_void_star &input_items,
                                              gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *)input_items[0];
      gr_complex *out = (gr_complex *)output_items[0];
      float *out2, *out3;
      if (output_items.size() >= 2)
        out2 = (float *)output_items[1];
      if (output_items.size() >= 3)
        out3 = (float *)output_items[2];
      int oidx = 0, iidx = 0;
      int ninp = ninput_items[0] - 3.0 * d_sps;
      if (ninp <= 0)
      {
        consume_each(0);
        return (0);
      }

      std::vector<tag_t> tags;
      get_tags_in_range(tags,
                        0,
                        nitems_read(0),
                        nitems_read(0) + ninp,
                        pmt::intern("time_est"));

      gr_complex sq,     // Squared input
          dly_conj,      // Input delayed sps and conjugated
          nlin_out,      // output of the nonlinearity
          in_interp;     // interpolated input
      float err_out = 0; // error output

      while (oidx < noutput_items && iidx < ninp)
      {
        // check to see if there's a tag to reset the timing estimate
        if (tags.size() > 0)
        {
          int offset = tags[0].offset - nitems_read(0);
          if ((offset >= iidx) && (offset < (iidx + d_sps)))
          {
            float center = (float)pmt::to_double(tags[0].value);
            if (center != center)
            { // test for NaN, it happens somehow
              tags.erase(tags.begin());
              goto out;
            }
            d_mu = center;
            iidx = offset;
            if (d_mu < 0)
            {
              d_mu++;
              iidx--;
            }
            d_div = 0;
            d_omega = d_sps;
            d_dly_conj_2 = d_dly_conj_1;
            // this keeps the block from outputting an odd number of
            // samples and throwing off downstream blocks which depend
            // on proper alignment -- for instance, a decimating FIR
            // filter.
            //                    if(d_div == 0 and d_osps == 2) oidx++;
            tags.erase(tags.begin());
          }
        }

      out:
        // the actual equation for the nonlinearity is as follows:
        // e(n) = in[n]^2 * in[n-sps].conj()^2
        // we then differentiate the error by subtracting the sample delayed by d_sps/2
        in_interp = d_interp->interpolate(&in[iidx], d_mu);
        sq = in_interp * in_interp;
        // conjugation is distributive.
        dly_conj = std::conj(d_dly_conj_2 * d_dly_conj_2);
        nlin_out = sq * dly_conj;
        // TODO: paper argues that some improvement can be had
        // if you either operate at >2sps or use a better numeric
        // differentiation method.
        err_out = std::real(nlin_out - d_dly_diff_1);
        if (d_div % 2)
        { // error loop calc once per symbol
          err_out = gr::branchless_clip(err_out, 3.0);
          d_omega += d_gain_omega * err_out;
          d_omega = d_sps + gr::branchless_clip(d_omega - d_sps, d_limit);
          d_mu += d_gain * err_out;
        }
        // output every other d_sps by default.
        if (!(d_div % 2) or d_osps == 2)
        {
          out[oidx] = in_interp;
          if (output_items.size() >= 2)
            out2[oidx] = err_out;
          if (output_items.size() >= 3)
            out3[oidx] = d_mu;
          oidx++;
        }
        d_div++;

        d_dly_conj_1 = in_interp;
        d_dly_conj_2 = d_dly_conj_1;
        d_dly_diff_1 = nlin_out;

        // update interpolator twice per symbol
        d_mu += d_omega;
        iidx += (int)floor(d_mu);
        d_mu -= floor(d_mu);
      }

      consume_each(iidx);
      return oidx;
    }

  } /* namespace ais */
} /* namespace gr */
