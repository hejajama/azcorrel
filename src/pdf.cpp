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

void PDF::PlotPdf(double Q)
{
    cout << "# PDF at Q=" << Q << "GeV" << endl;
    cout <<"# x     up    down   s  gluon" << endl;
    for (double x=1e-4; x<1; x*=1.1)
    {
        cout << x << " " << xq(x,Q,U) << " " <<  xq(x,Q,D) << " "
        << xq(x,Q,S) << " " << xq(x,Q,G) << endl;
    }

}
