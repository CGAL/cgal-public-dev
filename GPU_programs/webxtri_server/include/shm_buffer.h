// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/webxalci/include/shm_buffer.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:11 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef SHM_BUFFER_H
#define SHM_BUFFER_H

#include <qiodevice.h>
#include <qstring.h>

//! implements writing directly into shared memory region

class Shm_buffer : public QIODevice
{
public:
    Shm_buffer();
    Shm_buffer(char *, uint);
    ~Shm_buffer() { }
	
	bool setBuffer(char *, uint);

    bool  open( int );
    void  close();
    void  flush();

    Offset size() const;
    Offset at() const;
    bool  at( Offset );

    Q_LONG	  readBlock( char *p, Q_ULONG );
    Q_LONG	  writeBlock( const char *p, Q_ULONG );
    Q_LONG	  writeBlock( const QByteArray& data );
	      //{ return QIODevice::writeBlock(data); }
    Q_LONG	  readLine( char *p, Q_ULONG );

    int	  getch();
    int	  putch( int );
    int	  ungetch( int );

protected:
    //QByteArray a;
	char *shm_addr;
	uint max_size;
	
private:
    uint  a_size;
    
private:	// Disabled copy constructor and operator=
#if defined(Q_DISABLE_COPY)
    //#error "disabled copy constructor ?"
	Shm_buffer( const Shm_buffer & );
    Shm_buffer &operator=( const Shm_buffer & );
#endif
};

inline QIODevice::Offset Shm_buffer::size() const
{ return (Offset)a_size; }

inline QIODevice::Offset Shm_buffer::at() const
{ return ioIndex; }


#endif // SHM_BUFFER_H
