//Based on "root_v6.06.06.source.tar.gz":root-6.06.06/hist/hist/src/TH2.cxx, line 1330-1495 
//Double_t TH2D::KolmogorovTest(const TH1 *h2, Option_t *option) const

Double_t MyKolmogorovTest2D(const TH1 *h1, const TH1 *h2, const std::vector<int> vec1, const std::vector<int> vec2, Option_t *option)
{
   TString opt = option;
   opt.ToUpper();

   if (vec1.size()!=vec2.size()) {
     std::cout << "number of empty bin for x and y is inconsistent." << std::endl;
     return 1;
   }

   Double_t prb = 0;
   //TH1 *h1 = (TH1*)this;
   if (h2 == 0) return 0;
   const TAxis *xaxis1 = h1->GetXaxis();
   const TAxis *xaxis2 = h2->GetXaxis();
   const TAxis *yaxis1 = h1->GetYaxis();
   const TAxis *yaxis2 = h2->GetYaxis();
   Int_t ncx1 = xaxis1->GetNbins();
   Int_t ncx2 = xaxis2->GetNbins();
   Int_t ncy1 = yaxis1->GetNbins();
   Int_t ncy2 = yaxis2->GetNbins();

   // Check consistency of dimensions
   if (h1->GetDimension() != 2 || h2->GetDimension() != 2) {
      Error("KolmogorovTest","Histograms must be 2-D\n");
      return 0;
   }

   // Check consistency in number of channels
   if (ncx1 != ncx2) {
      Error("KolmogorovTest","Number of channels in X is different, %d and %d\n",ncx1,ncx2);
      return 0;
   }
   if (ncy1 != ncy2) {
      Error("KolmogorovTest","Number of channels in Y is different, %d and %d\n",ncy1,ncy2);
      return 0;
   }

   // Check consistency in channel edges
   Bool_t afunc1 = kFALSE;
   Bool_t afunc2 = kFALSE;
   Double_t difprec = 1e-5;
   Double_t diff1 = TMath::Abs(xaxis1->GetXmin() - xaxis2->GetXmin());
   Double_t diff2 = TMath::Abs(xaxis1->GetXmax() - xaxis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      Error("KolmogorovTest","histograms with different binning along X");
      return 0;
   }
   diff1 = TMath::Abs(yaxis1->GetXmin() - yaxis2->GetXmin());
   diff2 = TMath::Abs(yaxis1->GetXmax() - yaxis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      Error("KolmogorovTest","histograms with different binning along Y");
      return 0;
   }

   // Should we include Uflows, Oflows?
   Int_t ibeg = 1, jbeg = 1;
   Int_t iend = ncx1, jend = ncy1;
   if (opt.Contains("U")) {ibeg = 0; jbeg = 0;}
   if (opt.Contains("O")) {iend = ncx1+1; jend = ncy1+1;}

   Int_t i,j;
   Double_t sum1 = 0;
   Double_t sum2 = 0;
   Double_t w1   = 0;
   Double_t w2   = 0;
   for (i = ibeg; i <= iend; i++) {
      for (j = jbeg; j <= jend; j++) {

        bool skip_bin_flag = false;
        for(int iEmptyBin=0;iEmptyBin<(int)vec1.size();iEmptyBin++){
          if(i==vec1[iEmptyBin] && j==vec2[iEmptyBin]){
            skip_bin_flag = true;
            break;
          }
        }
        if(skip_bin_flag==true) continue;

         sum1 += h1->GetBinContent(i,j);
         sum2 += h2->GetBinContent(i,j);
         Double_t ew1 = h1->GetBinError(i,j);
         Double_t ew2 = h2->GetBinError(i,j);
         w1 += ew1*ew1;
         w2 += ew2*ew2;
         if(0){
           double tmp_entry = h1->GetBinContent(i,j);
           std::cout << "Bin(" << i << ", " << j << ") ==> Content = " << tmp_entry << ", pow(Content, 0.5) = " << pow(tmp_entry, 0.5) << std::endl
                     << "Error = " << h1->GetBinError(i,j) << std::endl;
	 }
      }
   }

   // Check that both scatterplots contain events
   if (sum1 == 0) {
      Error("KolmogorovTest","Integral is zero for h1=%s\n",h1->GetName());
      return 0;
   }
   if (sum2 == 0) {
      Error("KolmogorovTest","Integral is zero for h2=%s\n",h2->GetName());
      return 0;
   }
   // calculate the effective entries.
   // the case when errors are zero (w1 == 0 or w2 ==0) are equivalent to
   // compare to a function. In that case the rescaling is done only on sqrt(esum2) or sqrt(esum1)
   Double_t esum1 = 0, esum2 = 0;
   if (w1 > 0)
      esum1 = sum1 * sum1 / w1;
   else
      afunc1 = kTRUE;// use later for calculating z

   if (w2 > 0)
      esum2 = sum2 * sum2 / w2;
   else
      afunc2 = kTRUE;// use later for calculating z

   if (afunc2 && afunc1) {
      Error("KolmogorovTest","Errors are zero for both histograms\n");
      return 0;
   }

   // Find first Kolmogorov distance
   Double_t s1 = 1/sum1;
   Double_t s2 = 1/sum2;
   Double_t dfmax1 = 0;
   Double_t rsum1 = 0, rsum2 = 0;
   for (i=ibeg;i<=iend;i++) {
      for (j=jbeg;j<=jend;j++) {

        bool skip_bin_flag = false;
        for(int iEmptyBin=0;iEmptyBin<(int)vec1.size();iEmptyBin++){
          if(i==vec1[iEmptyBin] && j==vec2[iEmptyBin]){
            skip_bin_flag = true;
            break;
          }
        }
        if(skip_bin_flag==true) continue;

         rsum1 += s1*h1->GetBinContent(i,j);
         rsum2 += s2*h2->GetBinContent(i,j);
         dfmax1 = TMath::Max(dfmax1, TMath::Abs(rsum1-rsum2));
      }
   }

   // Find second Kolmogorov distance
   Double_t dfmax2 = 0;
   rsum1 = 0, rsum2 = 0;
   for (j=jbeg;j<=jend;j++) {
      for (i=ibeg;i<=iend;i++) {

        bool skip_bin_flag = false;
        for(int iEmptyBin=0;iEmptyBin<(int)vec1.size();iEmptyBin++){
          if(i==vec1[iEmptyBin] && j==vec2[iEmptyBin]){
            skip_bin_flag = true;
            break;
          }
        }
        if(skip_bin_flag==true) continue;

         rsum1 += s1*h1->GetBinContent(i,j);
         rsum2 += s2*h2->GetBinContent(i,j);
         dfmax2 = TMath::Max(dfmax2, TMath::Abs(rsum1-rsum2));
      }
   }

   // Get Kolmogorov probability: use effective entries, esum1 or esum2,  for normalizing it
   Double_t factnm;
   if (afunc1)      factnm = TMath::Sqrt(esum2);
   else if (afunc2) factnm = TMath::Sqrt(esum1);
   else             factnm = TMath::Sqrt(esum1*sum2/(esum1+esum2));

   // take average of the two distances
   Double_t dfmax = 0.5*(dfmax1+dfmax2);
   Double_t z = dfmax*factnm;
   if(0) std::cout << "z = " << z << std::endl;

   prb = TMath::KolmogorovProb(z);

   Double_t prb1 = 0, prb2 = 0;
   // option N to combine normalization makes sense if both afunc1 and afunc2 are false
   if (opt.Contains("N")  && !(afunc1 || afunc2 ) ) {
      // Combine probabilities for shape and normalization
      prb1   = prb;
      Double_t d12    = esum1-esum2;
      Double_t chi2   = d12*d12/(esum1+esum2);
      prb2   = TMath::Prob(chi2,1);
      // see Eadie et al., section 11.6.2
      if (prb > 0 && prb2 > 0) prb = prb*prb2*(1-TMath::Log(prb*prb2));
      else                     prb = 0;
   }

   // debug printout
   if (opt.Contains("D")) {
      printf(" Kolmo Prob  h1 = %s, sum1=%g\n",h1->GetName(),sum1);
      printf(" Kolmo Prob  h2 = %s, sum2=%g\n",h2->GetName(),sum2);
      printf(" Kolmo Probabil = %f, Max Dist = %g\n",prb,dfmax);
      if (opt.Contains("N"))
         printf(" Kolmo Probabil = %f for shape alone, =%f for normalisation alone\n",prb1,prb2);
   }
   // This numerical error condition should never occur:
   if (TMath::Abs(rsum1-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h1=%s\n",h1->GetName());
   if (TMath::Abs(rsum2-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h2=%s\n",h2->GetName());

   if(opt.Contains("M"))      return dfmax;// return avergae of max distance

   return prb;
}
