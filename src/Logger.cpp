//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


/** \file src/Logger.cpp
** \brief Message logging class.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#define __IN_LOGGER_CPP
#include<Logger.h>

// The global logging object
Logger Log;

Logger::Logger(void) :
    pDebugStream(&std::cout),
    pInfoStream(&std::cout),
    pErrorStream(&std::cerr),
    NullStream(NULL),
    Debugging(false)
{
    std::ios_base::sync_with_stdio(true);
}

void Logger::setupStream(std::ostream& Stream)
{
    Stream << std::setprecision(8) << std::setw(9) << std::left;
}

// Debug
void Logger::setDebugStream(std::ostream& Stream){ pDebugStream = &Stream; setupStream(Stream);}
std::ostream& Logger::debug(void){ return Debugging ? (*pDebugStream) : NullStream;}
std::ostream& debug(void){ return Log.debug(); }

void Logger::setDebugging(bool Debugging)
{
    this->Debugging = Debugging;
}

// Info
void Logger::setInfoStream(std::ostream& Stream){ pInfoStream = &Stream; setupStream(Stream);}
std::ostream& Logger::info(void){ return *pInfoStream;}
std::ostream& info(void){ return Log.info(); }

// Error
void Logger::setErrorStream(std::ostream& Stream){ pErrorStream = &Stream; setupStream(Stream);}
std::ostream& Logger::error(void){ return *pErrorStream;}
std::ostream& error(void){ return Log.error(); }
