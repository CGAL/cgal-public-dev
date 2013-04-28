// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.f
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
// File          : demos/webxalci/server.C
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:04 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#define NDEBUG 1

#include <signal.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/shm.h>
#include <pthread.h>
#include <semaphore.h>

#include <openssl/md5.h>
#include <time.h>

#include <arpa/inet.h>
#include <netdb.h>

#include <qapplication.h>
#include <qtimer.h>
#include <qpixmap.h>
#include <qpainter.h>
#include <qimage.h>
#include <qfile.h>

#include <fstream>
#include <list>

#define CGAL_NO_LEDA

#include <CGAL/Timer.h>
#include "include/shm_buffer.h"
#include "include/IPC.h"

#include "include/CGAL_includes.h"
#include "include/rasterizer.h"
#include "include/server.h"

//! external symbols declaration
extern XAlci_server *pxalci_server;
extern Rasterizer *prasterizer;
extern pthread_attr_t detach_attr;
extern pthread_mutex_t arcs_cache_mtx;
extern pthread_mutex_t time_mtx;
extern pthread_mutex_t comments_mtx;
extern pthread_mutex_t painter_mtx;
extern pthread_mutex_t active_cl_mtx;
extern pthread_mutex_t GPU_arithm_mtx;

extern pthread_mutex_t RS_alloc_mtx;
extern pthread_mutex_t NTL_arithm_mtx;
extern pthread_mutex_t CGAL_static_mtx;

extern pthread_mutex_t log_mtx;
extern pthread_cond_t active_job_cv;
extern sem_t ipc_msg_sem;
extern sem_t shadow_sem;

extern const char *rasterize_modes[];

//! this thread refreshes active task list and cancels timed out requests;
//! it also updates add comment delays
void *XAlci_server::multiplexer_thread(void *)
{

CGAL_precondition(switch off precondiotions!!);

    timespec ts;
    int i;
    while(1) {
    clock_gettime(CLOCK_REALTIME, &ts);
    ts.tv_sec += THREAD_LIST_REFRESH; 
    sem_timedwait(&shadow_sem, &ts);
    
//std::cout << "# of comments elements: " << comment_delay_list.size() << "\n";
    std::list<in_addr_t> to_remove_ip; // a list of entries to be removed
    pthread_mutex_lock(&comments_mtx);
    Comment_delay_list::const_iterator comment_it = comment_delay_list.begin();
    while(comment_it != comment_delay_list.end()) {
        if(comment_it->second >= ts.tv_sec) {
            comment_it++;
            continue;
        }
        to_remove_ip.push_back(comment_it->first);
        comment_it++;
    }
    std::list<in_addr_t>::iterator rm_ip_it = to_remove_ip.begin();
    while(rm_ip_it != to_remove_ip.end())
        comment_delay_list.erase(*rm_ip_it++);
    pthread_mutex_unlock(&comments_mtx);
        
    std::list<pthread_t> to_remove; // a list of entries to be removed
    pthread_mutex_lock(&time_mtx);
//std::cout << "refreshing thread_list: " << thread_list.size() << std::endl;
    Thread_list::iterator it = thread_list.begin();
    while(it != thread_list.end()) {
        if(it->second.timeout >= ts.tv_sec) {
            it++;
            continue;
        }
        Message_type mtype = it->second.msg_type;

        cancelled_id = it->first;
        if(mtype == ANALYSE) {
            write_log("Cancelling analysis request: PID = %d", it->first);
            pthread_mutex_lock(&arcs_cache_mtx);
            memcpy(&md_variable_set, &(it->second.md), sizeof(MD5_digest));
            active_job.erase(it->second.md);
            
        } else if(mtype == RASTERIZE) {
            write_log("Cancelling rasterize request: PID = %d", it->first);
        }

      // lock CORE internal mutexes otherwise the computation thread may be
      // cancelled while CORE mutex is still locked
        pthread_cancel(it->first);
        thread_cleanup_handler(&it->second);
        
        write_log("cancelling done");
        
        if(mtype == ANALYSE) {
        // notify all waiting threads that active job was removed, this will 
        // force them to exit with an error
            pthread_cond_broadcast(&active_job_cv);
            pthread_mutex_unlock(&arcs_cache_mtx);
        } 
        to_remove.push_back(it->first);
        it++;
    }
    
    std::list<pthread_t>::iterator it2 = to_remove.begin();
    while(it2 != to_remove.end())
        thread_list.erase(*it2++);
    pthread_mutex_unlock(&time_mtx);
    }
    return NULL;
}

typedef void *(*ALLOC_FUN) (size_t);
typedef void *(*REALLOC_FUN) (void *, size_t, size_t);
typedef void (*FREE_FUN) (void *, size_t);

ALLOC_FUN def_gmp_alloc, rs_push_alloc;
REALLOC_FUN def_gmp_realloc, rs_push_realloc;
FREE_FUN def_gmp_free, rs_push_free;

static pthread_t rs_heap_owner = -1u;

void* my_allocate(size_t s){

    if(rs_heap_owner != pthread_self()) {
        return def_gmp_alloc(s);
    }
    void *p;
    p = rs_push_alloc(s);

    return p;
}
 
void* my_reallocate(void *a, size_t o, size_t n){

//     printf("realloc: %x: %d - ", a, n);
    if(rs_heap_owner != pthread_self()) {
        return def_gmp_realloc(a, o, n);
    }
    void *p;
    p = rs_push_realloc(a, o, n);
    return p;
}

void my_free(void *a, size_t s){

    if(rs_heap_owner != pthread_self()) {
        def_gmp_free(a, s);
        return;
    }
    rs_push_free(a,s);
}

void my_push_func(ALLOC_FUN alloc_func,
    REALLOC_FUN realloc_func, FREE_FUN free_func) {

//     printf("My push func: %x %x %x\n",
//             alloc_func, realloc_func, free_func);

    if(rs_heap_owner != -1u) {
        printf("FATAL: attempt to acquire Rs heap while it is "
            "used by another thread: %d %d\n", rs_heap_owner,
                pthread_self());
        exit(1);
    }
    rs_push_alloc = alloc_func;
    rs_push_realloc = realloc_func;
    rs_push_free = free_func; 

//     pthread_mutex_lock(&CCPA_cache_mtx);
    rs_heap_owner = pthread_self();
//     printf("RS owner: %d\n", rs_heap_owner);
//     mp_set_memory_functions(my_allocate, my_reallocate, my_free);
}

void my_pop_func(void) {
//     printf("RS owner released\n", rs_heap_owner);
    rs_heap_owner = -1u;
//     mp_set_memory_functions(PTR_alloc, PTR_realloc, PTR_free);
//     pthread_mutex_unlock(&CCPA_cache_mtx);
}

void XAlci_server::thread_cleanup_handler(Thread_cleanup_info *info)
{
    // this indicates that an error occurred
    ((int *)info->shm_addr)[0] = 0xdeadbeef;
    ((int *)info->shm_addr)[1] = 0xdeadbeef;
    info->pmsg->shm_size = 8;
    
    if(info->host_addr != (in_addr_t)-1) { 
        pthread_mutex_lock(&active_cl_mtx);
        active_clients.erase(info->host_addr); 
        pthread_mutex_unlock(&active_cl_mtx);
    }
    
    if(msgsnd(mq_id, info->pmsg, sizeof(IPC_Message)-4, IPC_NOWAIT)== -1)
        err_msg("msgsnd");
    if(shmdt(info->shm_addr) == -1) // detach shared memory address
        err_msg("shmdt");
    write_log("cleanup handler done..");
}

Error_code XAlci_server::handle_analyse_request(Analysis_entry& analysis,
    char *ascii_poly, Thread_cleanup_info *info) {

    static int first_time = 0;
    
    Error_code ret = ERR_OK;
    Poly_int_vector poly_vec;
    if(ascii_poly != NULL) { // ASCII poly is nonzero meaning that
                             // MD5 checksum is not computed
        if(!parse_polynomial(ascii_poly, analysis.curve_ID, poly_vec))
            return ERR_INVALID_POLYNOMIAL;
    }
        
    pthread_mutex_lock(&arcs_cache_mtx);
    boost::multi_index::nth_index<Analysis_cache,1>::type& idx =
            analysis_cache.get<1>();
    boost::multi_index::nth_index_iterator<Analysis_cache,1>::type
            it = idx.find(analysis.curve_ID);
//     Analysis_cache::iterator it = analysis_cache.find(*it);

    pthread_mutex_lock(&comments_mtx);
    if(first_time == 0) { //! make sure this code executes only once

        mp_get_memory_functions(&def_gmp_alloc, &def_gmp_realloc,
                &def_gmp_free);

//         rs3_set_gmp_usr_memory_fncts(my_allocate, my_reallocate, my_free);
        mp_set_memory_functions(my_allocate, my_reallocate, my_free);

        set_push_gmp_alloc_fnct(my_push_func);
        set_pop_gmp_alloc_fnct(my_pop_func);
    }
    first_time = 1;
    pthread_mutex_unlock(&comments_mtx);
            
    if(it == idx.end()) {
        if(ascii_poly == NULL) {// identifier not found in the cache
            pthread_mutex_unlock(&arcs_cache_mtx);
            return ERR_RASTERIZE_GENERIC;
        }
        
        if(active_job.find(analysis.curve_ID) == active_job.end()) {
            std::cout << "new active job added\n";
            active_job.insert(analysis.curve_ID);
            pthread_mutex_unlock(&arcs_cache_mtx);

            // call actual curve analysis routine
            thread_list_insert(ANALYSE, &analysis.curve_ID, info);

            CGAL::Timer tm_anal;
            Points_2 points;
            try {
                tm_anal.start();
                prasterizer->analyse_curve(analysis, poly_vec);
                tm_anal.stop();

                std::cout << "######## analysis time: " << tm_anal.time() <<
                        "\n";
            }
            catch(...) {
                // catch whatever can be caught and report error
                if(cancelled_id == pthread_self())
                    throw; // rethrow exception if asynchronous cancel
                
                // catch whatever can be caught and report error
                thread_list_remove();
                write_log("internal server error\n");
                
                pthread_mutex_lock(&arcs_cache_mtx);
                active_job.erase(analysis.curve_ID);
                memcpy(&md_variable_set, &analysis.curve_ID,
                    sizeof(MD5_digest));
                pthread_cond_broadcast(&active_job_cv);
                pthread_mutex_unlock(&arcs_cache_mtx);
                return ERR_INVALID_DATA;
            }
            thread_list_remove();

            prasterizer->print_out(analysis);
            //timer.stop();
            pthread_mutex_lock(&arcs_cache_mtx);
            if(analysis_cache.size() >= MAX_ANALYSIS_CACHE_SIZE) {
                analysis_cache.pop_back();
                std::cout << "removing LRU entry from the cache\n";
            }
            analysis.one_curve = (poly_vec.size() == 1);
            analysis_cache.push_front(analysis);
            write_log("new entry added to the cache: %d",
                analysis_cache.size());
            memcpy(&md_variable_set, &analysis.curve_ID, sizeof(MD5_digest));
            active_job.erase(analysis.curve_ID);
            pthread_cond_broadcast(&active_job_cv);
        
        } else {
            write_log("waiting on condition variable");
            while(analysis.curve_ID != md_variable_set)
                pthread_cond_wait(&active_job_cv, &arcs_cache_mtx);
            // wait for condition variable
            idx = analysis_cache.get<1>();
            it = idx.find(analysis.curve_ID);
            if(it != idx.end()) {
                analysis = *it;
            } else {
                write_log("ERROR: analysis not found");
                ret = ERR_INVALID_DATA;
            }
        }
        
    } else { // mark this element as most recently used
        analysis = *it;
        analysis_cache.relocate(analysis_cache.project<0>(it),
                analysis_cache.begin());
    }
    pthread_mutex_unlock(&arcs_cache_mtx);
    //std::cout << "\nanalyse2 mutex UNlocked: " << n_locks << std::endl;
    return ret;
}

uint XAlci_server::handle_rasterize_request(const Analysis_entry& analysis,
    uint *indices, char *shm_addr, Thread_cleanup_info *info) {
    
    QPixmap plot(width, height);
    plot.fill();
    
    SHM_Data *data = (SHM_Data *)shm_addr;
    // copy field-by-field to avoid errors if fields are disordered
    CGAL::Bbox_2 box(data->x_min, data->y_min, data->x_max, data->y_max);

    std::string outp;
    // calling actual curve rendering routine 
    thread_list_insert(RASTERIZE, 0, info);
    switch(data->mode) {
    
    case DRAW_SUBDIV_1: {
        MD5_digest poly_hash;
        Poly_int_vector poly_vec;
        // parsing + testing for squarefreeness could be lengthy operation =>
        // must be controlled by timeouts
        if(!parse_polynomial((char *)indices, poly_hash, poly_vec)) {
            thread_list_remove();
            return -2u; // indicates a parser error
        }
        prasterizer->plot_subdivision(poly_vec, box, &plot);
        break;
    }        
    case LOCATE_POINT: {
        Point_coords *pcoords = (Point_coords *)((char *)indices +
            sizeof(MD5_digest));
        //std::cerr << pcoords->x << "; and " << pcoords->y << "\n";

        SHM_Point_query_reply *reply =
            (SHM_Point_query_reply *)shm_addr;

        bool res = prasterizer->locate_point(analysis, box, *pcoords, reply);
        thread_list_remove();
        return (res ? sizeof(SHM_Point_query_reply) : -1);
    }
    case HIGHLIGHT_VERTS:
        prasterizer->plot_points(analysis, indices, data->n_indices, box,
                data->mode, &plot);
        break;
        
    case HIGHLIGHT_FACES:
        prasterizer->plot_faces(analysis, indices, data->n_indices, box,
                data->mode, &plot);
        break;

    case HIGHLIGHT_ARCS:
    case DRAW_DEFAULT:
    case DRAW_FLAT:
        prasterizer->plot_arcs(analysis, indices, data->n_indices, box,
              data->mode, &plot);
        break;
                
    case PRINT_COORDS:
        prasterizer->plot_arcs(analysis, indices, data->n_indices, box,
              data->mode, &outp);

        break;

    case DRAW_GRAPH:
        prasterizer->draw_topology_graph(analysis, box, &plot);
        break;

    default:
        thread_list_remove();
        std::cerr << "bogus rendering mode: " << data->mode << std::endl;
        return -1;
    }
    thread_list_remove();        

    if(data->mode != PRINT_COORDS) {
        pthread_mutex_lock(&painter_mtx);
        Shm_buffer shm_buf(shm_addr, MAX_SEG_SIZE);
        shm_buf.open(IO_WriteOnly);
        plot.save(&shm_buf, "PNG");
        pthread_mutex_unlock(&painter_mtx);
        write_log("%d byte(s) written", shm_buf.size());
        return shm_buf.size();    
    } 
    
    strncpy(shm_addr, outp.c_str(), MAX_COORDS_PRINTOUT);
    //std::cerr << shm_addr << std::endl;
    int len = strlen(shm_addr);
    write_log("%d byte(s) written", len);
    return len;
}

void *XAlci_server::main_request_thread(IPC_Message *pmsg)
{
    IPC_Message ipc_msg;    
    int shm_id;
    char *shm_addr;
    memcpy(&ipc_msg, pmsg, sizeof(IPC_Message));
    // unlock semaphore, signalling that ipc_msg was successfully copied to
    // local memory
    sem_post(&ipc_msg_sem);
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
    
    key_t key = ipc_msg.shm_key;
    if((shm_id = shmget(key, MAX_SEG_SIZE, 0660)) == -1) {
        err_msg("shmget");
        return NULL;
    }
    //std::cout << "shared memory ID recevied.. " << shm_id << std::endl;
    if((int&)(shm_addr = (char *)shmat(shm_id, NULL, 0)) == -1) {
        err_msg("shmat");
        return NULL;
    }
        
    long type = ipc_msg.m_type, poly_size = 0;
    SHM_Data *shm_data = (SHM_Data *)shm_addr;
    MD5_digest poly_hash;
    uint n_indices = shm_data->n_indices, 
        *indices = (uint *)(shm_addr+sizeof(SHM_Data));
    pid_t cPID = shm_data->PID;
    in_addr host_addr = {shm_data->host_addr};
    char *poly = NULL;
    const char *pstr;
    bool invalid = false;
    
    // retrieve the time of request and a hostname of the client
    Request_info request_info;
    time_t time_result = time(NULL);
    pthread_mutex_lock(&time_mtx); // do we need this ?
    strncpy(request_info.time, asctime(localtime(&time_result)), 63);
    request_info.time[strlen(request_info.time)-1]='\0';
    // copy host ip-address since inet_ntoa returns a pointer to static data
    pstr = inet_ntoa(host_addr);
    strncpy(request_info.ip_address, pstr, 31);
    request_info.ip_address[31] = '\0';
    // copy the hostname
//     pstr = (const char *)gethostbyaddr(&host_addr, sizeof(in_addr), AF_INET);
//     pstr = (pstr == NULL ? (char *)"unrecognized host" :
//         ((hostent *)pstr)->h_name);
    pstr = "";
    strncpy(request_info.hostname, pstr, 255);
    request_info.hostname[255] = '\0'; 
    pthread_mutex_unlock(&time_mtx);

    switch(type) {
    case ANALYSE:
        ipc_msg.m_type = (cPID << 4)|ANALYSE_ACK;
        pstr = "CURVE_ANALYSE";
        poly = (char *)indices;
        poly_size = (int)ipc_msg.shm_size - sizeof(SHM_Data);
        break;
        
    case RASTERIZE:
        ipc_msg.m_type = (cPID << 4)|RASTERIZE_ACK;
        if(shm_data->mode == DRAW_SUBDIV_1) {
            pstr = "DRAW_SUBDIV_1";
            poly = (char *)indices;
            poly_size = (int)ipc_msg.shm_size - sizeof(SHM_Data);
        } else {
            memcpy(&poly_hash, (char *)indices + n_indices*sizeof(int),
                sizeof(MD5_digest));

            if(shm_data->mode >= DRAW_DEFAULT &&
                    shm_data->mode <= DRAW_GRAPH)
                pstr = rasterize_modes[shm_data->mode];
            else {
                pstr = (char *)"BOGUS_REQUEST";
                invalid = true;
            }
        }
        break;
        
    case COMMENT:
        ipc_msg.m_type = (cPID << 4)|COMMENT_ACK;
        pstr = "COMMENT";
        poly = (char *)indices;
        poly_size = (int)ipc_msg.shm_size - sizeof(SHM_Data);
        poly[poly_size] = '\0';
        break;
        
    default:
        pstr = "UNKNOWN";
        invalid = true;
    }
    write_log("%s: %s from %s (%s); pid: %d",
        request_info.time, pstr, request_info.ip_address,
            request_info.hostname, cPID);
    if(invalid)
        return NULL;
                    
    ipc_msg.err_code = ERR_TIMEOUT; // set timeout error code since computation
                                    // thread may be cancelled
    // install cleanup handler
    Thread_cleanup_info cleanup_info;
    cleanup_info.pmsg = &ipc_msg;
    cleanup_info.shm_addr = (void *)shm_addr;
    cleanup_info.host_addr = (in_addr_t)-1; // temporarily set to error value
//     pthread_cleanup_push(thread_cleanup_handler_proxy, &cleanup_info);
    unsigned do_cleanup = 1;
    
    if(type != COMMENT) {
        
        // wrong # of indices or wrong shm size
        if(n_indices > 2048||ipc_msg.shm_size > MAX_SEG_SIZE||poly_size < 0) { 
            write_log("Wrong client data..");
            ipc_msg.err_code = ERR_INVALID_DATA;
            goto Lexit;
        }
            
        pthread_mutex_lock(&active_cl_mtx);
        if(active_clients.size() > MAX_PENDING_REQUESTS) {
            pthread_mutex_unlock(&active_cl_mtx);
            write_log("Server overload; # of clients: %d",
                active_clients.size());
            ipc_msg.err_code = ERR_SERVER_OVERLOAD;
            goto Lexit;
        }
            
        boost::multi_index::nth_index_iterator<Active_clients,0>::type 
            it = active_clients.find(host_addr.s_addr);
        if(it != active_clients.end()) { // this client already has a request
/** *************************************************
**** HACK HACK HACK ******************************
****************************************************/
//             pthread_mutex_unlock(&active_cl_mtx);
//             write_log("Request pending..");
//             ipc_msg.err_code = ERR_REQUEST_PENDING;
//             goto Lexit;
/** *************************************************
**** HACK HACK HACK ******************************
****************************************************/

        } else { // otherwise make this client active
            active_clients.insert(Active_client_info(host_addr.s_addr));
            // from this point on we should keep track client's IP address
            cleanup_info.host_addr = host_addr.s_addr; 
        }
        pthread_mutex_unlock(&active_cl_mtx);
    
        if(poly != NULL)
            poly[poly_size]='\0';

        Analysis_entry analysis;
        if(type == ANALYSE || (type == RASTERIZE &&
                shm_data->mode != DRAW_SUBDIV_1)) {
            if(type == RASTERIZE)
                memcpy(&analysis.curve_ID, &poly_hash, sizeof(MD5_digest));
                
            ipc_msg.err_code = handle_analyse_request(analysis, poly,
                                            &cleanup_info);
            if(ipc_msg.err_code != ERR_OK) {
                goto Lexit;
            } 
        }

        // for ANALYSE request we send MD5 hash of a polynomial back to client
        if(type == ANALYSE) {

            SHM_Analysis_reply *reply = (SHM_Analysis_reply *)shm_addr;
            memcpy(&reply->curve_ID, &analysis.curve_ID, sizeof(MD5_digest));

            reply->n_faces = analysis.arr.number_of_faces();
            reply->n_edges = analysis.arr.number_of_edges();
            reply->n_vertices = analysis.arr.number_of_vertices();
            reply->n_isolated = analysis.arr.number_of_isolated_vertices();
            /////// ATTENTION!! the size now is one byte more ?

            strcpy(&reply->print_out, analysis.arcs_print.c_str());
            reply->vidx_start = analysis.arcs_print.length() + 1; // 0-char

            strcpy((&reply->print_out) + reply->vidx_start,
                analysis.verts_print.c_str());
            reply->fidx_start = reply->vidx_start +
                analysis.verts_print.length() + 1;

            strcpy((&reply->print_out) + reply->fidx_start,
                analysis.faces_print.c_str());

            ipc_msg.shm_size = sizeof(SHM_Analysis_reply) +
                reply->fidx_start + analysis.faces_print.length()+1;
          
        } else if(type == RASTERIZE) {
            if(shm_data->mode != DRAW_SUBDIV_1 &&
                    shm_data->mode != DRAW_GRAPH)
                write_log("# of indices: %d", n_indices);

            // in case of SUBDIV_1 poly = indices (zero character is already
            // inserted): poly[poly_size] = '\0'; returns -2 in case of parser
            // error
            ipc_msg.shm_size = handle_rasterize_request(analysis, indices,
                shm_addr, &cleanup_info);
            if((int)ipc_msg.shm_size < 0) {
                ipc_msg.err_code = (ipc_msg.shm_size == -1 ?
                    ERR_RASTERIZE_GENERIC : ERR_INVALID_POLYNOMIAL);
                goto Lexit;
            }
        }
        
        // remove client from the active list
        pthread_mutex_lock(&active_cl_mtx);
        active_clients.erase(host_addr.s_addr); 
        pthread_mutex_unlock(&active_cl_mtx);
        
    } else if(!handle_comment_request(&request_info, poly, host_addr)) {
        ipc_msg.err_code = ERR_REQUEST_PENDING;
        goto Lexit;
    }

    do_cleanup = 0; // do not execute cleanup handler
    ipc_msg.err_code = ERR_OK; // no errors
    if(msgsnd(mq_id, &ipc_msg, sizeof(ipc_msg)-4, IPC_NOWAIT) == -1)
        perror("msgsnd");
    if(shmdt(shm_addr) == -1) 
        err_msg("shmdt");

Lexit:
    if(do_cleanup) {
        thread_cleanup_handler(&cleanup_info);
    }
//     pthread_cleanup_pop(do_cleanup); 
    return NULL;
}


void XAlci_server::run()
{
    // spawn multiplexor thread
    if(pthread_create(&mplex_id, &detach_attr, multiplexer_thread_proxy, 
            NULL) != 0) {
        err_msg("multiplexer pthread_create");
        return;
    }

#if CGAL_BISOLVE_USE_GPU_RESULTANTS
    CGAL::internal::GPU_gcd& obj1 = CGAL::internal::GPU_gcd::instance();
    CGAL::internal::GPU_resultant& obj2 =
             CGAL::internal::GPU_resultant::instance();

    (void)obj1;
    (void)obj2;

#endif

    while(1) {
        std::cout << server_id << ": waiting for clients..." << std::endl;
        IPC_Message ipc_msg;
        if(msgrcv(mq_id, &ipc_msg, sizeof(ipc_msg)-4, MSGS_SERVER, 0) == -1) {
            err_msg("msgrcv");
            err_exit();
        } 
        if(ipc_msg.m_type == ANALYSE||ipc_msg.m_type == RASTERIZE||
            ipc_msg.m_type == COMMENT) {
            pthread_t child;
            if(pthread_create(&child, &detach_attr, main_request_thread_proxy,
                    &ipc_msg)==0) {
                // decrement semaphore count to prevent simultaneous 
                // acccess to ipc_msg structure by another client when it is
                // still not copied to thread's local memory
                sem_wait(&ipc_msg_sem);
            } else
                perror("pthread_create");
            
        } else { 
            std::cout << "Unknown message type: " << ipc_msg.m_type 
                    << std::endl;
        }
    }
}

#define STDIN           0       /* what is stdin's file descriptor? */
#define STDOUT          1       /*   stdout's? */
#define STDERR          2       /*   stderr's? */

// edit /etc/init.d/xalci_server
//install using: update-rc.d xalci_server defaults

// NOTE NOTE NOTE: for debugging use: set follow-fork-mode child

int main(int argc, char *argv[]) 
{
    switch (fork()) {
    case -1:
        std::cout << "cannot fork..\n\n";
        exit(0);
        break;
    case 0:
        int fd;
        /* child process */
        if((fd = open("/dev/null", O_RDWR, 0)) >= 0) {
            (void) dup2(fd, STDIN);
//           (void) dup2(fd, STDOUT);
//           (void) dup2(fd, STDERR);
            if(fd != STDERR)
               (void) close(fd);
        }
        break;
    default:
        /* parent process should just die */
        exit(0);
    }

    int ret = 1, inst;
    for(inst = 1; inst < N_SERVER_INSTANCES; inst++) {
        ret = fork();
        if(ret == 0) // break for child process
            break;
    }

    if(ret == 0) {
        std::cout << "hello I am child process: " << getpid()
            << "; inst: " << inst << "\n";
//         exit(0);
    } else
        inst = 0;
    std::cout << "Parent process: " << getpid() << " continue.. " <<
         inst << "\n";

    umask(0);

    /* Create a new SID for the child process */
    int sid = setsid();
    if(sid < 0) {
        /* Log any failure */
        std::cout << "cannot setsid..\n\n";
        exit(1);
    }

    static XAlci_server server_obj(inst);   
    static Rasterizer rasterizer_obj;
    QApplication q_app(argc, argv); 
    pxalci_server = &server_obj;
    prasterizer = &rasterizer_obj;
    
    pxalci_server->setup(700, 525);
    pxalci_server->run();
    return 0;
}
