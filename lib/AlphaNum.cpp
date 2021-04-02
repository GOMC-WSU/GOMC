#include "AlphaNum.h"

AlphaNum::AlphaNum(){}

/* Alpha numberic A, B, ... Z, AA, AB, ... AZ, */
/* Index from 0 to A */

std::string AlphaNum::uint2String(uint stringSuffix){

    std::stringstream ss;
    char charSuffix;
    int intermediate = stringSuffix;
    int remainder = 0;
    do {
        charSuffix = 'A';
        remainder = intermediate % 26;
        /* Increment Char A until reach suffix or 27 which will be Z. */
        for (int j = 0; j < remainder; j++){
            charSuffix++;
        }
        ss << charSuffix;
        intermediate /= 26;
        intermediate--;
    } while (intermediate >= 0);

    std::string backwards = ss.str();
    std::string forwards;

    return forwards.assign(backwards.rbegin(), backwards.rend());
}