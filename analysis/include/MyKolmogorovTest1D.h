//Based on "root_v6.06.06.source.tar.gz":root-6.06.06/hist/hist/src/TH1.cxx, line 7609-7765
//Double_t TH1::KolmogorovTest(const TH1 *h2, Option_t *option) const

Double_t MyKolmogorovTest1D(const TH1 *h1, const TH1 *h2, const std::vector<int> vec1, Option_t *option)
{
   TString opt = option;
   opt.ToUpper();

   Double_t prob = 0;
   //TH1 *h1 = (TH1*)this;
   if (h2 == 0) return 0;
   const TAxis *axis1 = h1->GetXaxis();
   const TAxis *axis2 = h2->GetXaxis();
   Int_t ncx1   = axis1->GetNbins();
   Int_t ncx2   = axis2->GetNbins();

   // Check consistency of dimensions
   if (h1->GetDimension() != 1 || h2->GetDimension() != 1) {
      Error("KolmogorovTest","Histograms must be 1-D\n");
      return 0;
   }

   // Check consistency in number of channels
   if (ncx1 != ncx2) {
      Error("KolmogorovTest","Number of channels is different, %d and %d\n",ncx1,ncx2);
      return 0;
   }

   // empty the buffer. Probably we could add as an unbinned test
   //if (fBuffer) ((TH1*)this)->BufferEmpty();

   // Check consistency in channel edges
   Double_t difprec = 1e-5;
   Double_t diff1 = TMath::Abs(axis1->GetXmin() - axis2->GetXmin());
   Double_t diff2 = TMath::Abs(axis1->GetXmax() - axis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      Error("KolmogorovTest","histograms with different binning");
      return 0;
   }

   Bool_t afunc1 = kFALSE;
   Bool_t afunc2 = kFALSE;
   Double_t sum1 = 0, sum2 = 0;
   Double_t ew1, ew2, w1 = 0, w2 = 0;
   Int_t bin;
   Int_t ifirst = 1;
   Int_t ilast = ncx1;
   // integral of all bins (use underflow/overflow if option)
   if (opt.Contains("U")) ifirst = 0;
   if (opt.Contains("O")) ilast = ncx1 +1;
   for (bin = ifirst; bin <= ilast; bin++) {

     bool skip_bin_flag = false;
      for(int iEmptyBin=0;iEmptyBin<(int)vec1.size();iEmptyBin++){
        if(bin==vec1[iEmptyBin]){
          skip_bin_flag = true;
          break;
        }
      }
      if(skip_bin_flag==true) continue;

      sum1 += h1->GetBinContent(bin);
      sum2 += h2->GetBinContent(bin);
      ew1   = h1->GetBinError(bin);
      ew2   = h2->GetBinError(bin);
      w1   += ew1*ew1;
      w2   += ew2*ew2;
      if(0){
         double tmp_entry = h1->GetBinContent(bin);
         std::cout << "Bin " << bin << " ==> Content = " << tmp_entry << ", pow(Content, 0.5) = " << pow(tmp_entry, 0.5) << std::endl
                   << "Error = " << h1->GetBinError(bin) << std::endl;
       }
   }
   if (sum1 == 0) {
      Error("KolmogorovTest","Histogram1 %s integral is zero\n",h1->GetName());
      return 0;
   }
   if (sum2 == 0) {
      Error("KolmogorovTest","Histogram2 %s integral is zero\n",h2->GetName());
      return 0;
   }

   // calculate the effective entries.
   // the case when errors are zero (w1 == 0 or w2 ==0) are equivalent to
   // compare to a function. In that case the rescaling is done only on sqrt(esum2) or sqrt(esum1)
   Double_t esum1 = 0, esum2 = 0;
   if (w1 > 0)
      esum1 = sum1 * sum1 / w1;
   else
      afunc1 = kTRUE;    // use later for calculating z

   if (w2 > 0)
      esum2 = sum2 * sum2 / w2;
   else
      afunc2 = kTRUE;    // use later for calculating z

   if (afunc2 && afunc1) {
      Error("KolmogorovTest","Errors are zero for both histograms\n");
      return 0;
   }


   Double_t s1 = 1/sum1;
   Double_t s2 = 1/sum2;

   // Find largest difference for Kolmogorov Test
   Double_t dfmax =0, rsum1 = 0, rsum2 = 0;

   for (bin=ifirst;bin<=ilast;bin++) {
      rsum1 += s1*h1->GetBinContent(bin);
      rsum2 += s2*h2->GetBinContent(bin);
      dfmax = TMath::Max(dfmax,TMath::Abs(rsum1-rsum2));
   }

   // Get Kolmogorov probability
   Double_t z, prb1=0, prb2=0, prb3=0;

   // case h1 is exact (has zero errors)
  if  (afunc1)
      z = dfmax*TMath::Sqrt(esum2);
  // case h2 has zero errors
  else if (afunc2)
      z = dfmax*TMath::Sqrt(esum1);
  else
     // for comparison between two data sets
     z = dfmax*TMath::Sqrt(esum1*esum2/(esum1+esum2));

   prob = TMath::KolmogorovProb(z);

   // option N to combine normalization makes sense if both afunc1 and afunc2 are false
   if (opt.Contains("N") && !(afunc1 || afunc2 ) ) {
      // Combine probabilities for shape and normalization,
      prb1 = prob;
      Double_t d12    = esum1-esum2;
      Double_t chi2   = d12*d12/(esum1+esum2);
      prb2 = TMath::Prob(chi2,1);
      // see Eadie et al., section 11.6.2
      if (prob > 0 && prb2 > 0) prob *= prb2*(1-TMath::Log(prob*prb2));
      else                      prob = 0;
   }
   //// X option. Pseudo-experiments post-processor to determine KS probability
   //const Int_t nEXPT = 1000;
   //if (opt.Contains("X") && !(afunc1 || afunc2 ) ) {
   //   Double_t dSEXPT;
   //   TH1 *hExpt = (TH1*)(gDirectory ? gDirectory->CloneObject(this,kFALSE) : gROOT->CloneObject(this,kFALSE));
   //   // make nEXPT experiments (this should be a parameter)
   //   prb3 = 0;
   //   for (Int_t i=0; i < nEXPT; i++) {
   //      hExpt->Reset();
   //      hExpt->FillRandom(h1,(Int_t)esum2);
   //      dSEXPT = KolmogorovTest(hExpt,"M");
   //      if (dSEXPT>dfmax) prb3 += 1.0;
   //   }
   //   prb3 /= (Double_t)nEXPT;
   //   delete hExpt;
   //}

   // debug printout
   if (opt.Contains("D")) {
      printf(" Kolmo Prob  h1 = %s, sum bin content =%g  effective entries =%g\n",h1->GetName(),sum1,esum1);
      printf(" Kolmo Prob  h2 = %s, sum bin content =%g  effective entries =%g\n",h2->GetName(),sum2,esum2);
      printf(" Kolmo Prob     = %g, Max Dist = %g\n",prob,dfmax);
      if (opt.Contains("N"))
         printf(" Kolmo Prob     = %f for shape alone, =%f for normalisation alone\n",prb1,prb2);
      //if (opt.Contains("X"))
      //   printf(" Kolmo Prob     = %f with %d pseudo-experiments\n",prb3,nEXPT);
   }
   // This numerical error condition should never occur:
   if (TMath::Abs(rsum1-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h1=%s\n",h1->GetName());
   if (TMath::Abs(rsum2-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h2=%s\n",h2->GetName());

   if(opt.Contains("M"))      return dfmax;
   else if(opt.Contains("X")) return prb3;
   else                       return prob;
}
