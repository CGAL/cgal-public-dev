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
// File          : demos/webxalci/shm_buffer.C
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:04 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================


#include "include/shm_buffer.h"
//#include "include/xalci.h"

//#include <stdlib.h>
//#include "qfile.h"

Shm_buffer::Shm_buffer()
{
    setFlags( IO_Direct );
    shm_addr = NULL;
	max_size = 0;
	ioIndex = 0;
}

Shm_buffer::Shm_buffer(char *shm_addr_, uint max_size_) : 
		shm_addr(shm_addr_), max_size(max_size_)
{
    setFlags( IO_Direct );
    a_size = 0;
    ioIndex = 0;
}

bool Shm_buffer::setBuffer(char *shm_addr_, uint max_size_) 
{
    if ( isOpen() ) {
        return FALSE;
    }
    shm_addr = shm_addr_;
    a_size = 0;
	max_size = max_size_;
	ioIndex = 0;
    return TRUE;
}

bool Shm_buffer::open( int m  )
{
    if ( isOpen() ) {                           // buffer already open
        return FALSE;
    }
    setMode( m );
    if ( m & IO_Truncate ) {                    // truncate buffer
        a_size = 0;
    }
    if ( m & IO_Append ) {                      // append to end of buffer
        ioIndex = a_size;
    } else {
        ioIndex = 0;
    }
    setState( IO_Open );
    resetStatus();
    return TRUE;
}

void Shm_buffer::close()
{
    if ( isOpen() ) {
        setFlags( IO_Direct );
        ioIndex = 0;
    }
}

void Shm_buffer::flush()
{
    return;
}

bool Shm_buffer::at( Offset pos )
{
#if defined(QT_CHECK_STATE)
    if ( !isOpen() ) {
        qWarning( "Shm_buffer::at: Buffer is not open" );
        return FALSE;
    }
#endif
    if ( pos > a_size ) {
        return FALSE;
    }
    ioIndex = pos;
    return TRUE;
}

Q_LONG Shm_buffer::readBlock( char *p, Q_ULONG len )
{
#if defined(QT_CHECK_STATE)
    if ( !p ) {
	qWarning( "Shm_buffer::readBlock: Null pointer error" );
	return -1;
    }
    if ( !isOpen() ) {                          // buffer not open
        qWarning( "Shm_buffer::readBlock: Buffer not open" );
        return -1;
    }
    if ( !isReadable() ) {                      // reading not permitted
        qWarning( "Shm_buffer::readBlock: Read operation not permitted" );
        return -1;
    }
#endif
    if ( ioIndex + len > a_size ) {   // overflow
        if(ioIndex >= a_size) 
            return 0;
        else 
            len = a_size - ioIndex;
    }
    memcpy(p, shm_addr + ioIndex, len);
    ioIndex += len;
    return len;
}

Q_LONG Shm_buffer::writeBlock( const char *p, Q_ULONG len )
{
    if ( len == 0 )
        return 0;

#if defined(QT_CHECK_NULL)
    if ( p == 0 ) {
        qWarning( "Shm_buffer::writeBlock: Null pointer error" );
        return -1;
    }
#endif
#if defined(QT_CHECK_STATE)
    if ( !isOpen() ) {                          // buffer not open
        qWarning( "Shm_buffer::writeBlock: Buffer not open" );
        return -1;
    }
    if ( !isWritable() ) {                      // writing not permitted
        qWarning( "Shm_buffer::writeBlock: Write operation not permitted" );
        return -1;
    }
#endif
    if(ioIndex + len > max_size) {              // overflow
        setStatus( IO_ResourceError );
        return -1;
    }
    memcpy(shm_addr+ioIndex, p, len);
    ioIndex += len;
    if(a_size < ioIndex)
		a_size = ioIndex;
	return len;
}


/*!
  \reimp
*/

Q_LONG Shm_buffer::readLine( char *p, Q_ULONG maxlen )
{
#if defined(QT_CHECK_NULL)
    if ( p == 0 ) {
        qWarning( "Shm_buffer::readLine: Null pointer error" );
        return -1;
    }
#endif
#if defined(QT_CHECK_STATE)
    if ( !isOpen() ) {                          // buffer not open
        qWarning( "Shm_buffer::readLine: Buffer not open" );
        return -1;
    }
    if ( !isReadable() ) {                      // reading not permitted
        qWarning( "Shm_buffer::readLine: Read operation not permitted" );
        return -1;
    }
#endif
    if ( maxlen == 0 )
        return 0;
    Q_ULONG start = ioIndex;
    char *d = shm_addr + ioIndex;
    maxlen--;                                   // make room for 0-terminator
    if(a_size - ioIndex < maxlen)
        maxlen = a_size - ioIndex;
    while ( maxlen-- ) {
        if ( (*p++ = *d++) == '\n' )
            break;
    }
    *p = '\0';
    ioIndex = d - shm_addr;
    return ioIndex - start;
}

int Shm_buffer::getch()
{
#if defined(QT_CHECK_STATE)
    if ( !isOpen() ) {                          // buffer not open
        qWarning( "Shm_buffer::getch: Buffer not open" );
        return -1;
    }
    if ( !isReadable() ) {                      // reading not permitted
        qWarning( "Shm_buffer::getch: Read operation not permitted" );
        return -1;
    }
#endif
    if ( ioIndex+1 > a_size ) {               // overflow
        setStatus( IO_ReadError );
        return -1;
    }
    return uchar(*(shm_addr+ioIndex++));
}

int Shm_buffer::putch( int ch )
{
#if defined(QT_CHECK_STATE)
    if ( !isOpen() ) {                          // buffer not open
        qWarning( "Shm_buffer::putch: Buffer not open" );
        return -1;
    }
    if ( !isWritable() ) {                      // writing not permitted
        qWarning( "Shm_buffer::putch: Write operation not permitted" );
        return -1;
    }
#endif
    if(ioIndex + 1 > max_size) {                // overflow
        char buf[1];
        buf[0] = (char)ch;
        if ( writeBlock(buf,1) != 1 )
            return -1;                          // write error
    } else {
        *(shm_addr + ioIndex++) = (char)ch;
        if( a_size < ioIndex )
            a_size = ioIndex;
    }
    return ch;
}

int Shm_buffer::ungetch( int ch )
{
#if defined(QT_CHECK_STATE)
    if ( !isOpen() ) {                          // buffer not open
        qWarning( "Shm_buffer::ungetch: Buffer not open" );
        return -1;
    }
    if ( !isReadable() ) {                      // reading not permitted
        qWarning( "Shm_buffer::ungetch: Read operation not permitted" );
        return -1;
    }
#endif
    if ( ch != -1 ) {
        if ( ioIndex )
            ioIndex--;
        else
            ch = -1;
    }
    return ch;
}



