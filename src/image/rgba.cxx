// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: rgba.cxx,v 1.1.1.1 2003-02-12 16:51:52 whitaker Exp $

#include "image/rgba.h"

rgba operator*(const rgba& rgba_in,  byte value )
{
    rgba new_rgba;
    new_rgba._color[R] = rgba_in._color[R]*value;
    new_rgba._color[G] = rgba_in._color[G]*value;
    new_rgba._color[B] = rgba_in._color[B]*value;
    new_rgba._color[A] = rgba_in._color[A]*value;
    return(new_rgba);
}

rgba operator*( byte value, const rgba& rgba_in)
{
    rgba new_rgba;
    new_rgba._color[R] = rgba_in._color[R]*value;
    new_rgba._color[G] = rgba_in._color[G]*value;
    new_rgba._color[B] = rgba_in._color[B]*value;
    new_rgba._color[A] = rgba_in._color[A]*value;
    return(new_rgba);

}

rgba operator*(const rgba& rgba_in, float value)
{
    rgba new_rgba;
    new_rgba._color[R] = (byte)(rgba_in._color[R]*value + 0.5);
    new_rgba._color[G] = (byte)(rgba_in._color[G]*value + 0.5);
    new_rgba._color[B] = (byte)(rgba_in._color[B]*value + 0.5);
    new_rgba._color[A] = (byte)(rgba_in._color[A]*value + 0.5);
    return(new_rgba);

}

rgba operator*(float value, const rgba& rgba_in)
{
    rgba new_rgba;
    new_rgba._color[R] = (byte)(rgba_in._color[R]*value + 0.5);
    new_rgba._color[G] = (byte)(rgba_in._color[G]*value + 0.5);
    new_rgba._color[B] = (byte)(rgba_in._color[B]*value + 0.5);
    new_rgba._color[A] = (byte)(rgba_in._color[A]*value + 0.5);
    return(new_rgba);
}

