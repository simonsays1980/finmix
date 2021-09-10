/******************************************************************************
*
* Copyright (C) 2013 Lars Simon Zehnder. All Rights Reserved.
*
* Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
*
* This file is part of the R package 'finmix'.
*
* 'finmix' is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published
* by the Free Software Foundatio, either version 3 of the License, or
* any later version.
*
* 'finmix' is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with 'finmix'. If not, see <http://www.gnu.org/licenses/>.
*
******************************************************************************/


#ifndef BASE_H
#define BASE_H

#include <RcppArmadillo.h>

// ==============================================================
// BASE class (to be reviewed)
// --------------------------------------------------------------
/**
 * @brief   Base class for all mixin layers.
 * @detail  This is the base class for all mixin layers defined.
 *          In particular it defines next to constructor and
 *          destructor the virtual member functions 'BASE::update'
 *          and 'BASE::store' to be implemented or redefined by
 *          inheriting classes.
 * @see FIX, IND, POST, HIER, ADAPTER
 * @author Lars Simon Zehnder
 *
 * ==============================================================
 * @review  As any combination of layers begins with the FIX
 *          the latter one has no super class and needs thereby
 *          also no base class, as it defines 'Node' and
 *          'Output' classes as well as 'update()' and 'store()'
 *          methods. The BASE class is probably only needed, if
 *          there is no real hierarchy to the system of layers.
 *          The same seems to hold for the adapter.
 * --------------------------------------------------------------
 */
class BASE {
public:
BASE ()
{
}
virtual ~BASE ()
{
}
/*
 * Function to update parameters.
 * Specified in all classes inheriting
 * from BASE
 */
virtual void update()
{
}
/*
 * Function to store values. Specified
 * all classes inheriting from BASE
 */
virtual void store(const unsigned int&)
{
}
};
#endif
