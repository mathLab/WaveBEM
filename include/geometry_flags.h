//---------------------------------------------------------------------------
//    $Id: fe_update_flags.h,v 1.31 2005/10/24 04:33:03 guido Exp $
//    Version: $Name:  $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__geometry_flags_h
#define __deal2__geometry_flags_h

#include <deal.II/base/config.h>
/**
 * The enum type given to the constructors of LocalAssembleBase objects,
 * telling those objects which data to assemble on each mesh cell.
 * When the GlobalAssembler calls the local one, it checks for each flag,
 * and if it finds one, it assemble the corresponding object.
 *
 * By default, all flags are off, i.e. no procedure will be called.
 *
 * You can select more than one flag by concatenation
 * using the bitwise or operator|(GeometryFlags,GeometryFlags).
 */
enum GeometryFlags
{
  //! No update
  none                      = 0,
  // Surface nodes
  water         = 0x0001,
  boat          = 0x0002,
  walls         = 0x0004,
  inflow          = 0x0008,
  // Edge nodes
  near_water          = 0x0010,
  near_boat         = 0x0020,
  near_walls          = 0x0040,
  edge          = 0x0080,
  keel          = 0x0100,
  near_inflow         = 0x0200,
  // Position characterization
  right_side                = 0x0400,
  left_side                 = 0x0800,
  // transom on water or boat
  transom_on_boat           = 0x1000,
  transom_on_water          = 0x2000,
  // pressure condition nodes
  pressure                  = 0x4000,
  near_pressure                  = 0x8000
};




/**
 * Global operator which returns an object in which all bits are set
 * which are either set in the first or the second argument. This
 * operator exists since if it did not then the result of the bit-or
 * <tt>operator |</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * GeometryFlags.
 */
inline
GeometryFlags
operator | (GeometryFlags f1, GeometryFlags f2)
{
  return static_cast<GeometryFlags> (
           static_cast<unsigned int> (f1) |
           static_cast<unsigned int> (f2));
}

/**
 * Global operator which returns an object in which all bits are set
 * which are either set in the first or the second argument. This
 * operator exists since if it did not then the result of the bit-or
 * <tt>operator |</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * GeometryFlags.
 */
inline
GeometryFlags
operator ^ (GeometryFlags f1, GeometryFlags f2)
{
  return static_cast<GeometryFlags> (
           static_cast<unsigned int> (f1) ^
           static_cast<unsigned int> (f2));
}

/**
 * Global operator which sets the bits from the second argument also
 * in the first one.
 */
inline
GeometryFlags &
operator |= (GeometryFlags &f1, GeometryFlags f2)
{
  f1 = f1 | f2;
  return f1;
}

/**
 * Global operator which sets the bits from the second argument also
 * in the first one only if they are not both true.
 */
inline
GeometryFlags &
operator ^= (GeometryFlags &f1, GeometryFlags f2)
{
  f1 = f1 ^ f2;
  return f1;
}


/**
 * Global operator which returns an object in which all bits are set
 * which are set in the first as well as the second argument. This
 * operator exists since if it did not then the result of the bit-and
 * <tt>operator &</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * GeometryFlags.
 */
inline
GeometryFlags
operator & (GeometryFlags f1, GeometryFlags f2)
{
  return static_cast<GeometryFlags> (
           static_cast<unsigned int> (f1) &
           static_cast<unsigned int> (f2));
}


/**
 * Global operator which clears all the bits in the first argument if
 * they are not also set in the second argument.
 */
inline
GeometryFlags &
operator &= (GeometryFlags &f1, GeometryFlags f2)
{
  f1 = f1 & f2;
  return f1;
}



#endif
