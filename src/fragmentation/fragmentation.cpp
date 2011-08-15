/*
 * Virtual class to hide different fragmentation functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <string>
#include "../config.hpp"
#include "fragmentation.hpp"

FragmentationFunction::FragmentationFunction()
{


}

std::string FragmentationFunction::GetString()
{
    return "not specified";
}
