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

#ifndef XALCI_CLIENT_H
#define XALCI_CLIENT_H

//! unique key in order to randomize IPC message queue ID
#define WEBXALCI_UNIQUE_KEY  0xdeadbeef
// a filename to be used to generate a unique IPC key
#define KEY_FILENAME        "/etc/security/limits.conf" 
    //"/root/xalci/problem_curve.txt"
   
// CAUTION!!! this parameter must be changed very carefully
// taking into account maximal size of client-server messages
#define MAX_SEG_SIZE        256*1024L    // max shared segment size

#define CLIENT_TIMEOUT 280 // timeout for client messages
#define ADD_COMMENT_DELAY 60

// server name used to check CGI-referrer
#define SERVER_NAME     "http://localhost:8080"

//! number of server copies to fork at the beginning
#define N_SERVER_INSTANCES 1

// maximum and minimum window dimensions
#define WINDOW_DIM_MAX  1e10
#define WINDOW_DIM_MIN  1e-13

#define BUF_SIZE 40*1024


/**
ANALYSE request: xalci.cgi?id=1&curve=bla
RASTERIZE request: xalci.cgi?id=2&all=<whether to draw all arcs>&
    r=<Rasterize_mode>&curveid=<curveid>&xmin=...&ymax=...
COMMENT request: xalci.cgi?id=3&...

*/

//! message type: (PID << 4) | MSG_TYPE - to identify a process
enum Message_type {
    ANALYSE = 1,       // analyse curve
    RASTERIZE = 2,     // rasterize curve
    COMMENT = 3,       // send a comment
    MSGS_SERVER = -3,  // all messages dedicated for server
    ANALYSE_ACK = 4,   // analyse request server response
    RASTERIZE_ACK = 5, // rasterize server response
    COMMENT_ACK = 6    // send comment server response
};

//! rasterization mode
enum Rasterize_mode {
    BOGUS = -1,
    DRAW_DEFAULT = 0,       // segment renderer (default)
    DRAW_FLAT = 1,          // rasterize complete curve in black color
    DRAW_SUBDIV_1 = 2,      // use 1D subdivision (experimental)
    HIGHLIGHT_ARCS = 3,     // highlight arcs
    HIGHLIGHT_VERTS = 4,    // highlight vertices
    HIGHLIGHT_FACES = 5,    // highlight faces
    LOCATE_POINT = 6,       // locate point on the arrangement
    PRINT_COORDS = 7,       // prints out the point coordinates
    DRAW_GRAPH = 8          // draws a topology graph of a curve
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
                                  // requests processed >= MAX_CLIENTS)
    ERR_RASTERIZE_GENERIC = -10   // generic error during rasterize request
};

//! type of client requests
enum Request_type {
    REQ_DEFAULT = 0,         // no request: display web-page only
    REQ_ANALYSE = 1,         // analyse request
    REQ_RASTERIZE = 2,       // rasterize default
    REQ_COMMENT = 3          // send comment
    
    /*REQ_SUBDIV = 3,          // rasterize subdivision
    REQ_HIGHLIGHT_ARCS = 5,  // draw features highlighten
    REQ_HIGHLIGHT_VERTS = 6,  // draw features highlighten
    REQ_HIGHLIGHT_FACES = 7,  // draw features highlighten*/
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
    Rasterize_mode mode;    //! rasterization mode
    long n_indices;         //! number of indices (zero for analyse requests)
    double x_min;           //! dimensions of a drawing window (rasterize only)
    double x_max;
    double y_min;
    double y_max;
};

struct Point_coords {
    int x, y;
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

//! describes a reply format for curve analysis request
struct SHM_Analysis_reply {
    MD5_digest curve_ID; // hash identifying the set of polynomials
    int n_faces;         // # of faces in arrangement
    int n_edges;         // # of edges
    int n_vertices;      // # of vertices total
    int n_isolated;      // # of isolated vertices

// arcs print-out (starting from &print_out) followed by vertices print-out,
// followed by faces print-out
    int vidx_start;      // starting position of vertices print-out
    int fidx_start;      // starting position of faces print-out
    char print_out;      // marks the beginning of arcs print-out
};
// arcs: from &print_out; size=vidx_start
// verts: from print_out+vidx_start, sise=fidx_start-vidx_start
// faces from print_out+fidx_start till the end

enum {QUERY_ON_POINT, QUERY_ON_EDGE, QUERY_ON_FACE };

//! describes a reply format for point location query
struct SHM_Point_query_reply {
    int type;   // feature type
    int index;  // index within arrangement
};

//! global functions
void msg(const char *format, ...);
void echo(const char *text);
inline void err_msg(const char *text);
void err_exit();
void sigint_handler(int);

class CGI_Client 
{
public:

    CGI_Client() : x_min(-2.0), x_max(2.0), y_min(-1.5), y_max(1.5), 
            err_code(ERR_OK), shm_id(-1), shm_addr(NULL), arcs_list(NULL),
            verts_list(NULL), faces_list(NULL), mq_id(-1)
    {
        memset(&poly_hash, 0, sizeof(MD5_digest)); 
    }
    
    ~CGI_Client();
    
    bool process(Request_type req, Rasterize_mode mode, int server_id);
    void web_interface(Request_type req, int server_id);
    void run();
    bool get_params(Request_type req, Rasterize_mode mode = BOGUS);
    void comment_sent();
    bool show_dummy();
    
    std::string parse_output(char *poly);
       
    int get_msg_queue_id() const
    { return mq_id; }
    
protected:
    double x_min, x_max, y_min, y_max; // dimensions of drawing window
    Point_coords location;
    
    MD5_digest poly_hash;   // polynomial MD5 hash
    uint time_stamp;        // a time-stamp of the last access from cookies

    int n_faces, n_edges;   // arrangement data
    int n_vertices, n_isolated;

    Error_code err_code;       // error code from the server
    
    int shm_id;         // shared memory segment ID for IPC communications
    char *shm_addr;     // shared memory start address

    char *arcs_list, *verts_list, *faces_list;
    int arcs_size, verts_size, faces_size;

    int PID;
    int checkbox_all;   // stores the value of checkboxes
    
    int mq_id;          // message queue ID to communicate with a server
  
};

#endif

