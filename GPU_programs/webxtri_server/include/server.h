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
// File          : demos/webxalci/include/server.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:11 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*! \file server.h
 *  \brief Defines class \c XTri_server 
 *  
 *  Main server class: handling IPC requests from clients, sending commands to
 *  \c Skeletonizer. \file skeletonizer.h must be included *before* this header
 */

#ifndef XTRI_SERVER_H
#define XTRI_SERVER_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>

using boost::multi_index::multi_index_container;
using boost::multi_index::get;
using boost::multi_index::project;

//! maximal number of cache entries in curve analysis cache
#define MAX_ANALYSIS_CACHE_SIZE  500     
//! timeout for curve analysis phase (in seconds)
#define ANALYSIS_TIMEOUT    110
//! timeout for triangulation phase (in seconds)
#define TRIANGULATE_TIMEOUT   110
//! specifies a refresh time for active thread list (in seconds)
#define THREAD_LIST_REFRESH 4
//! maximal number of simultaneous pending requests allowed
#define MAX_PENDING_REQUESTS 15
//! anti-spam delay to send comments (in seconds)
#define ADD_COMMENT_DELAY 60

//! number of server copies to fork at the beginning
#define N_SERVER_INSTANCES 5

std::ostream& operator <<(std::ostream& out, const MD5_digest& md);

//! \brief main xalci server class
class XTri_server 
{
public:
    //!\name public typedefs
    //!@{

    //! hash function for \c Analysis_cache
    struct MD5_hash_func {
        size_t operator()(const MD5_digest& key) const {
            return (key.md[0]+key.md[1]+key.md[2]+key.md[3]);
        }
    };
    
    //! set of unique surface IDs
    typedef std::set< MD5_digest, MD5_compare > Surface_IDs;

    //! stores a set of host IPs whose requests are currently being processed
    //! by the server
    typedef boost::multi_index::multi_index_container<
        Active_client_info,
        boost::multi_index::indexed_by<
            boost::multi_index::hashed_unique<
            BOOST_MULTI_INDEX_MEMBER(Active_client_info,in_addr_t,host_addr)
    > > > 
    Active_clients;
    
    //! a list of computation threads processing client requests (required
    //! to control thread timeouts)
    typedef std::map<pthread_t, Thread_cleanup_info> Thread_list;

    //! stores a set of host IPs who has recently sent comments (to prevent
    //! server spamming)
    typedef std::map<in_addr_t, int> Comment_delay_list;

    //!@}
public:
    //!\name public methods
    //!@{

    XTri_server(int _id) : server_id(_id)
    { }

    //! server setup
    void setup();
    
    //! timer thread: contols request timeouts and cancels pending requests
    void *multiplexer_thread(void *);
    
    //! main client service thread: spawned each time the server gets a request
    void *main_request_thread(IPC_Message *pmsg);
    
    //! curve analysis request handling
    Error_code handle_analyse_request(MD5_digest& surface_ID,
         char *ascii_poly, Thread_cleanup_info *info);
    //! triangulation request handling
    uint handle_triangulate_request(const MD5_digest& surface_ID,
        SHM_Data *data, void *triangle_data, uint max_data_sz,
        Thread_cleanup_info *info);

    //! writes a comment to the file
    bool handle_comment_request(Request_info *request_info, char *body,
        in_addr host_addr);
    
    //! cleanup handler for client threads
    void thread_cleanup_handler(Thread_cleanup_info *info);
    
    //! inserts the calling thread id into the thread list indicating pending
    //! computations
    void thread_list_insert(Message_type type, MD5_digest *pmd,
                            Thread_cleanup_info *info);

    //! removes the calling thread from the thread list if a lengthy
    //! computation was successful
    void thread_list_remove();
    
    //! \brief parses polynomial in ASCII format, returns \c false in case of
    //! error
    //!
    //! in case of success returns precached squarefree part of the polynomial 
    bool parse_polynomial(char *ascii_poly, MD5_digest& poly_hash,
        Polynomial_3& poly);
            
    //! server listening loop
    void run();
               
    //! writes a log file 
    void write_log(const char *fmt, ...);
    
    //! destructor
    ~XTri_server();
    
    //!@}
private:
    //!\name private data members
    //!@{

    //! stores a list of curves currently being analyised
    std::set<MD5_digest, MD5_compare> active_job; 
    
    //! ??
    MD5_digest md_variable_set;
    
    int mq_id;          //! message queue ID to communicate with a server
    int server_id;      //! ID of this server instance 
    pthread_t mplex_id; //! ID of multiplexer thread
    pthread_t cancelled_id; //! ID of a thread being cancelled

    Surface_IDs surface_IDs; //! MD5 checksums of analysed surfaces
    //Analysis_cache analysis_cache; //! surface analysis cache
    
    //! stores a set of host IPs whose requests are currently being processed
    //! by the server
    Active_clients active_clients; 
    //! a list of computation threads processing client requests (required
    //! to control thread timeouts)
    Thread_list thread_list; 
    //! a list of IP-addresses who has recently sent add comment request 
    //! (to protect the server from spam)
    Comment_delay_list comment_delay_list; 
    //! log and comment files output streams
    std::ofstream logfile, comment_file;
    
    //!@}
};

#endif // XTRI_SERVER_H

