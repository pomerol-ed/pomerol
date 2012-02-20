//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2012 Igor Krivenko <igor@shg.ru>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


/** \file src/Logger.h
** \brief Message logging class.
**
** \author Igor Krivenko (igor@shg.ru)
*/
#ifndef __INCLUDE_LOGGER_H
#define __INCLUDE_LOGGER_H

#include<iostream>
#include<iomanip>

class Logger {
    std::ostream* pDebugStream;
    std::ostream* pInfoStream;
    std::ostream* pErrorStream;

    std::ostream NullStream;

    bool Debugging;

    void setupStream(std::ostream& Stream);

public:
    Logger(void);

    void setDebugStream(std::ostream& Stream);
    void setInfoStream(std::ostream& Stream);
    void setErrorStream(std::ostream& Stream);

    void setDebugging(bool Debugging);

    std::ostream& debug(void);
    std::ostream& info(void);
    std::ostream& error(void);
};

std::ostream& debug(void);
std::ostream& info(void);
std::ostream& error(void);

#define MSG_PREFIX            __FILE__ << ":" << __LINE__ << ": "
#define DEBUG(MSG)            debug() << MSG_PREFIX << MSG << std::endl;
#define INFO(MSG)             info() << MSG << std::endl;
#define INFO_NONEWLINE(MSG)   info() << MSG << std::flush;
#define ERROR(MSG)            error() << MSG_PREFIX << MSG << std::endl;

#ifndef __IN_LOGGER_CPP
extern Logger Log;
#endif // endif :: #ifdef __IN_LOGGER_CPP

#endif // endif :: #ifndef __INCLUDE_HDF5STORAGE_H
