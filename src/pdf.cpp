/*
 * Virtual class to hide different parton distribution functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "pdf.hpp"

void PDF::Initialize(int param)
{
}

std::string PDF::GetString()
{
    return "not specified";
}

PDF::~PDF()
{
}

void PDF::PlotPdf(REAL Q)
{
    cout << "# PDF at Q=" << Q << "GeV" << endl;
    cout <<"# x     up    down " << endl;
    for (REAL x=1e-4; x<1; x*=1.1)
    {
        cout << x << " " << xq(x,Q,U) << " " <<  xq(x,Q,D) << endl;
    }

}
