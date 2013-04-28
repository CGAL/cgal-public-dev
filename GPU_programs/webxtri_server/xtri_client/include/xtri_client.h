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
// Library       : 
// File          : xalci_client.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:11 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*!\file xalci_client.h
 *  CGI client for \c xalci_server application. Communication with server is
 *  done through linux IPC
 */

#ifndef XTRI_CLIENT_H
#define XTRI_CLIENT_H

//! unique key in order to randomize IPC message queue ID
#define WEBXTI_UNIQUE_KEY  0
// a filename to be used to generate a unique IPC key
#define KEY_FILENAME        "/etc/security/time.conf"
    //"/root/xalci/problem_curve.txt"
   
// CAUTION!!! this parameter must be changed very carefully
// taking into account maximal size of client-server messages
#define MAX_SEG_SIZE        512*1024L    // max shared segment size

#define CLIENT_TIMEOUT 240 // timeout for client messages
#define PING_TIMEOUT   2   // timeout for ping messages

#define ADD_COMMENT_DELAY 60

// server name used to check CGI-referrer
#define SERVER_NAME     "http://localhost:8080"

#define BUF_SIZE 40*1024

//! number of server copies to fork at the beginning
#define N_SERVER_INSTANCES 5

/**
ANALYSE request: xalci.cgi?id=1&curve=bla
RASTERIZE request: xalci.cgi?id=2&all=<whether to draw all arcs>&
    r=<Rasterize_mode>&curveid=<curveid>&xmin=...&ymax=...
COMMENT request: xalci.cgi?id=3&...

*/

//! message type: (PID << 4) | MSG_TYPE - to identify a process
enum Message_type {
    ANALYSE = 1,        // analyse surface
    TRIANGULATE = 2,    // compute triangulation
    COMMENT = 3,       // send a comment
    PING = 4,          // check for server connection
    MSGS_SERVER = -5,      // all messages dedicated for server
    ANALYSE_ACK = 4,   // analyse request server response
    TRIANGULATE_ACK = 5, // triangulate server response
    COMMENT_ACK = 6,    // send comment server response
    PING_ACK = 7        // ping response
};

//! type of client requests
enum Request_type {
    REQ_ANALYSE = 1,         // analyse request
    REQ_TRIANGULATE = 2,     // triangulate request
    REQ_COMMENT = 3,          // send comment
    REQ_PING = 4
};

//! triangulation mode
enum Triangulate_mode {
    BOGUS = -1,
    TRIANGULATE_DEFAULT = 0,       // default
    TRIANGULATE_ABS_BOUNDS = 1     // triangulate using absolute bounds
};

enum Error_code {
    ERR_OK = 0,                   // no errors
    ERR_INVALID_POLYNOMIAL = -1,  // invalid polynomial format
    ERR_INVALID_DATA = -2,        // invalid client parameters
    ERR_TIMEOUT = -3,             // analyse/rasterize phase timeout
    ERR_INVALID_REFERRER = -6,    // an attempt to access the script from
                                  // outside the server
    ERR_SERVER_TIMEOUT = -7,      // no connection to the server
    ERR_REQUEST_PENDING = -8,     // a request from this IP is already being
                                  // processed by the server
    ERR_SERVER_OVERLOAD = -9,     // the server is overloaded (number of
                                  // requests processed >= MAX_CLIENTS
    ERR_TRIANGULATE_GENERIC = -10   // generic error during rasterize request
};

//! describes the format of message queue
struct IPC_Message {   
    long m_type;             // request type
    union {
        key_t shm_key;       // the key of a shared memory region
        Error_code err_code; // error code: 0 indicates no errors
                             // (server messages only)
    };
    long shm_size;           // the size of a shared memory region
};

// shm format:
// MSG_ANALYSE:   PID (4), n_indices(4) == 0, < polynomial string >
// MSG_RASTERIZE: PID (4), n_indices(4) != 0, < array of indices >, 
//      < polynomial string >

//! describes the format of a shared memory region
struct SHM_Data {
    pid_t PID;              //! client's PID
    in_addr_t host_addr;    //! client's IP-address
    Triangulate_mode mode;    //! triangulation mode (obsolete ?)
    //! skeletonizer specific parameters
    double en_left, en_right;
    double en_bottom, en_top;
    double z_below, z_above;
    int sx, sy;
};

//! thread info to control timeout when the server is unreachable
struct Thread_info {
    IPC_Message *pmsg;  
    uint response;
};

//! MD5 checksum to uniquely identify polynomials
struct MD5_digest
{
    MD5_digest() { }
    bool operator ==(const MD5_digest& md) const
    {
        return (memcmp(this, &md, sizeof(MD5_digest))==0);
    }
    bool operator !=(const MD5_digest& md) const
    {
        return (!operator ==(md));
    }
    unsigned md[4];
};

//! describes a reply format for surface analysis request
struct SHM_Analysis_reply {
    MD5_digest surface_ID; // hash identifying the set of polynomials
};

//! global functions
void msg(const char *format, ...);
void echo(const char *text);
void err_msg(const char *text);
void err_exit();
void sigint_handler(int);

class CGI_Client 
{
public:

    CGI_Client() : err_code(ERR_OK), shm_id(-1), shm_addr(NULL), mq_id(-1)
    {
        memset(&poly_hash, 0, sizeof(MD5_digest)); 
    }
    
    ~CGI_Client();
    
    bool process(Request_type req, Triangulate_mode mode = BOGUS,
                        int server_id = 0);
    void run();
    bool get_params(Request_type req, Triangulate_mode mode = BOGUS);
    void comment_sent();
    bool show_dummy();

    void send_error_code();
    
    std::string parse_output(char *poly);
       
    int get_msg_queue_id() const
    { return mq_id; }
    
protected:
    MD5_digest poly_hash;   // polynomial MD5 hash
    uint time_stamp;        // a time-stamp of the last access from cookies

    Error_code err_code;       // error code from the server
    
    int shm_id;         // shared memory segment ID for IPC communications
    char *shm_addr;     // shared memory start address

    SHM_Data in_shm_data;
    int PID;
    int mq_id;          // message queue ID to communicate with a server
  
};

#endif

