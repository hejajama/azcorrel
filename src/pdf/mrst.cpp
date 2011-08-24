/*
 * Wrapper class to use mrst99 PDF
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "mrst.hpp"
#include "mrst99.h"
#include <string>

/*
 * Default value of param is -1
 */
void MRST::Initialize(int param)
{
    set=param;
    if (param==-1) set=1;
    mrst = new c_mrst();
}

MRST::~MRST()
{
    delete mrst;
}

double MRST::xq(double x, double q, Parton p)
{
    mrst->mrst99(x, q, set);
    switch(p)
    {
        case UVAL:
            return mrst->cont.upv;
        case DVAL:
            return mrst->cont.dnv;
        case USEA:
            return mrst->cont.usea;
        case DSEA:
            return mrst->cont.dsea;
        case G:
            return mrst->cont.glu;
        case U:
            return xq(x,q,UVAL) + xq(x,q,USEA);
        case D:
            return xq(x,q,DVAL) + xq(x,q,USEA);
        case S:
            return mrst->cont.str;
        case B:
            return mrst->cont.bot;
        case C:
            return mrst->cont.chm;
        default:
            cerr << "Parton " << p << " is not implemented! " << LINEINFO << endl;
    };
    return 0;
}

std::string MRST::GetString()
{
    return "MRST99";
}
