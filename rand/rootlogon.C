// rootlogon for ShUtil
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 63;
    // const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.80, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
    gStyle->SetNumberContours( NCont );

    gInterpreter->AddIncludePath("./ShUtil");
    gInterpreter->AddIncludePath("./macros");
    gInterpreter->AddIncludePath("./include");
    gInterpreter->AddIncludePath("./src");
    std::cout << "Adding the include path to : "  << gSystem->GetIncludePath( ) << std::endl;

}
