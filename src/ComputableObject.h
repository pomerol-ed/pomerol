/** \file src/ComputableObject.h
** \brief A parent abstract class for all complex computable objects, i.e. Greens Function, Density matrix, Two Particle GF etc.
** 
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#include "Misc.h"

#ifndef __INCLUDE_COMPUTABLEOBJECT_H
#define __INCLUDE_COMPUTABLEOBJECT_H

/** This abstract class is a prototype for any computable object in the current code, i.e. Greens Function, Two Particle GF, Density Matrix etc...
 *  It defines the current state of calculation and declares virtual methods "compute" and "dump". 
 *  The enumeration of statuses is listed in the Misc.h file
 */
class ComputableObject{
protected:
    unsigned short Status;                      //!< Current status of an object
public:
    unsigned short getStatus(){return Status;}; //!< Get the current status of an object   
    ComputableObject(){Status=Constructed;};              //!< Constructor
//    virtual void prepare(){};                   //!< Prepare all the containers of the object, do preliminary fast routines
//    virtual void compute(){};                   //!< Do the most expensive calculation. Finishing it means the object has finished all calculation jobs
};

#endif
