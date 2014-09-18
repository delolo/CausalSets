/*
 * StringTools.h
 *
 *  Created on: 24 Feb 2014
 *      Author: mb909
 */

#ifndef STRINGTOOLS_H_
#define STRINGTOOLS_H_

double SigDigs(double dValue, int nSignificantDigits) {
    double dSign = (dValue > 0.0) ? 1 : -1;
    dValue *= dSign;

    int nOffset = static_cast<int>(log10(dValue));
    if (nOffset > 0) ++nOffset;

    double dScale = pow(10.0, nSignificantDigits - nOffset);

    return dSign * static_cast<double>(
            round(dValue * dScale)
            ) / dScale;
}


#endif /* STRINGTOOLS_H_ */
