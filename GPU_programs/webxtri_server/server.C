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
#include <fstream>
#include <list>

#include "include/IPC.h"

#include "include/CGAL_includes.h"
#include "include/skeletonizer.h"
#include "include/server.h"

//! external symbols declaration
extern XTri_server *pxtri_server;
extern Skeletonizer *pskeletonizer;
extern pthread_attr_t detach_attr;
extern pthread_mutex_t SFA_cache_mtx;
extern pthread_mutex_t algorithms_mtx;
extern pthread_mutex_t time_mtx;
extern pthread_mutex_t comments_mtx;
extern pthread_mutex_t active_cl_mtx;
extern pthread_mutex_t GPU_arithm_mtx;

extern pthread_mutex_t RS_alloc_mtx;
extern pthread_mutex_t NTL_arithm_mtx;
extern pthread_mutex_t CGAL_static_mtx;

extern pthread_mutex_t log_mtx;
extern pthread_cond_t active_job_cv;
extern sem_t ipc_msg_sem;
extern sem_t shadow_sem;

//! this thread refreshes active task list and cancels timed out requests;
//! it also updates add comment delays
void *XTri_server::multiplexer_thread(void *)
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
            write_log("Cancelling analysis request: PID = %u", it->first);
            pthread_mutex_lock(&SFA_cache_mtx);
            memcpy(&md_variable_set, &(it->second.md), sizeof(MD5_digest));
            active_job.erase(it->second.md);
            
        } else if(mtype == TRIANGULATE) {
            write_log("Cancelling triangulate request: PID = %u", it->first);
        }

        pthread_cancel(it->first);
        // it's dirty solution but we cannot do better
        thread_cleanup_handler(&it->second);
        
        write_log("cancelling done");
        
        if(mtype == ANALYSE) {
        // notify all waiting threads that active job was removed, this will 
        // force them to exit with an error
            pthread_cond_broadcast(&active_job_cv);
            pthread_mutex_unlock(&SFA_cache_mtx);
        } 
        to_remove.push_back(it->first);
        it++;
    }

    // TODO: need to clear active threads list as well!!
    
    std::list<pthread_t>::iterator it2 = to_remove.begin();
    while(it2 != to_remove.end()) {
        write_log("removing PID: %u\n", *it2);

        thread_list.erase(*it2++);
    }
    pthread_mutex_unlock(&time_mtx);
    }
    return NULL;
}

void XTri_server::thread_cleanup_handler(Thread_cleanup_info *info)
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

Error_code XTri_server::handle_analyse_request(MD5_digest& surface_ID,
    char *ascii_poly, Thread_cleanup_info *info) {
    
    static int first_time = 0;
    
    Error_code ret = ERR_OK;
    Polynomial_3 input_poly;
    if(ascii_poly != NULL) { // ASCII poly is nonzero meaning that
                             // MD5 checksum is not computed
        if(!parse_polynomial(ascii_poly, surface_ID, input_poly))
            return ERR_INVALID_POLYNOMIAL;
    }

    pthread_mutex_lock(&SFA_cache_mtx);
    Surface_IDs::iterator it = surface_IDs.find(surface_ID);

    if(it == surface_IDs.end()) {
        if(ascii_poly == NULL) {// identifier not found in the cache
            pthread_mutex_unlock(&SFA_cache_mtx);
            return ERR_TRIANGULATE_GENERIC;
        }
        
        if(active_job.find(surface_ID) == active_job.end()) {
            std::cout << "new active job added\n";
            active_job.insert(surface_ID);
            pthread_mutex_unlock(&SFA_cache_mtx);

            // call actual curve analysis routine
            thread_list_insert(ANALYSE, &surface_ID, info);

            CGAL::Timer tm_anal;
            try {
                tm_anal.start();
                if(!pskeletonizer->analyse_surface(surface_ID, input_poly))
                    throw "exception"; // well not the best solution..
                tm_anal.stop();

                std::cout << "######## analysis time: " << tm_anal.time() <<
                        "\n";
            }
            catch(...) {
                pthread_mutex_unlock(&algorithms_mtx);
                
                // catch whatever can be caught and report error
                if(cancelled_id == pthread_self())
                    throw; // rethrow exception if asynchronous cancel
                
                thread_list_remove();
                write_log("internal server error\n");
                
                pthread_mutex_lock(&SFA_cache_mtx);
                active_job.erase(surface_ID);
                memcpy(&md_variable_set, &surface_ID, sizeof(MD5_digest));
                pthread_cond_broadcast(&active_job_cv);
                pthread_mutex_unlock(&SFA_cache_mtx);
                return ERR_INVALID_DATA;
            }
            thread_list_remove();

            //timer.stop();
            pthread_mutex_lock(&SFA_cache_mtx);
            surface_IDs.insert(surface_ID);
            
            write_log("new entry added to the cache: %d", surface_IDs.size());
            memcpy(&md_variable_set, &surface_ID, sizeof(MD5_digest));
            active_job.erase(surface_ID);
            pthread_cond_broadcast(&active_job_cv);
        
        } else {
            write_log("waiting on condition variable");
            while(surface_ID != md_variable_set)
                pthread_cond_wait(&active_job_cv, &SFA_cache_mtx);
            // wait for condition variable
            it = surface_IDs.find(surface_ID);
            if(it == surface_IDs.end()) {
                write_log("ERROR: analysis not found");
                ret = ERR_INVALID_DATA;
            }
        }
    }
    // make this surface current if it is already in the cache
    if(it != surface_IDs.end() && !pskeletonizer->load_surface(surface_ID)) {
        ret = ERR_TRIANGULATE_GENERIC;
    }
    pthread_mutex_unlock(&SFA_cache_mtx);
    return ret;
}

uint XTri_server::handle_triangulate_request(const MD5_digest& surface_ID,
        SHM_Data *data, void *triangle_data, uint max_data_sz,
        Thread_cleanup_info *info) {

    thread_list_insert(TRIANGULATE, 0, info);

    bool use_auto_bounds = true;
    
    switch(data->mode) {
    case TRIANGULATE_DEFAULT: {
        use_auto_bounds = true;
        break;
    }
    case TRIANGULATE_ABS_BOUNDS:
        use_auto_bounds = false;
        break;
        
    default:
        thread_list_remove();
        std::cerr << "bogus rendering mode: " << data->mode << std::endl;
        return -1;
    }
    uint sz = pskeletonizer->skeletonize(surface_ID, data, triangle_data,
             max_data_sz, use_auto_bounds);
    thread_list_remove();
    return sz;
}

void *XTri_server::main_request_thread(IPC_Message *pmsg)
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

    char *poly = 0;
    pid_t cPID = shm_data->PID;
    in_addr host_addr = {shm_data->host_addr};
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

    pstr = "";
    strncpy(request_info.hostname, pstr, 255);
    request_info.hostname[255] = '\0'; 
    pthread_mutex_unlock(&time_mtx);

    switch(type) {
    case ANALYSE:
        ipc_msg.m_type = (cPID << 4)|ANALYSE_ACK;
        pstr = "SURFACE_ANALYSE";
        poly = (char *)(shm_addr + sizeof(SHM_Data));
        poly_size = (int)ipc_msg.shm_size - sizeof(SHM_Data);
        break;
        
    case TRIANGULATE:
        ipc_msg.m_type = (cPID << 4)|TRIANGULATE_ACK;
        if(shm_data->mode != TRIANGULATE_DEFAULT) {
            pstr = (char *)"BOGUS_REQUEST";
            invalid = true;
        }
        pstr = "TRIANGULATE_DEFAULT";
        memcpy(&poly_hash, shm_addr + sizeof(SHM_Data),
                sizeof(MD5_digest));
        break;

    case COMMENT:
        ipc_msg.m_type = (cPID << 4)|COMMENT_ACK;
        pstr = "COMMENT";
        break;

    case PING:
        ipc_msg.m_type = (cPID << 4)|PING_ACK;
        pstr = "PING";
        break;
        
    default:
        pstr = "UNKNOWN";
        invalid = true;
    }
//     write_log("%s: %s from %s (%s); PID: %u",
//         request_info.time, pstr, request_info.ip_address,
//             request_info.hostname, cPID);
    write_log("%s: %s from %s; PID: %u",
         request_info.time, pstr, request_info.ip_address, cPID);
    if(invalid)
        return NULL;
                    
    ipc_msg.err_code = ERR_TIMEOUT; // set timeout error code since computation
                                    // thread may be cancelled
    // install cleanup handler
    Thread_cleanup_info cleanup_info;
    cleanup_info.pmsg = &ipc_msg;
    cleanup_info.shm_addr = (void *)shm_addr;
    cleanup_info.host_addr = (in_addr_t)-1; // temporarily set to error value
    //pthread_cleanup_push(thread_cleanup_handler_proxy, &cleanup_info);
    unsigned do_cleanup = 1;
    
    if(type != COMMENT && type != PING) {
        
        // wrong # of indices or wrong shm size
        if(ipc_msg.shm_size > MAX_SEG_SIZE || poly_size < 0) {
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
    
        if(poly != 0)
            poly[poly_size]='\0';

        MD5_digest surface_ID;
        if(type == TRIANGULATE)
            memcpy(&surface_ID, &poly_hash, sizeof(MD5_digest));
        
        ipc_msg.err_code = handle_analyse_request(surface_ID, poly,
                              &cleanup_info);
        if(ipc_msg.err_code != ERR_OK) {
            goto Lexit;
        }

        SHM_Analysis_reply *reply = (SHM_Analysis_reply *)shm_addr;
        void *triangle_data = (void *)(shm_addr + sizeof(SHM_Analysis_reply));

        uint max_data_sz, tri_data_sz;
        try {
            max_data_sz = MAX_SEG_SIZE - 16 - sizeof(SHM_Analysis_reply);
            tri_data_sz = handle_triangulate_request(surface_ID, shm_data,
                    triangle_data, max_data_sz, &cleanup_info);
        } catch(...) {
            pthread_mutex_unlock(&algorithms_mtx);
            if(cancelled_id == pthread_self())
                throw; // rethrow exception if asynchronous cancel

            thread_list_remove();
            write_log("triangulate: internal server error\n");
            ipc_msg.err_code = ERR_INVALID_DATA;
            goto Lexit;
        }

        if((int)tri_data_sz < 0) {
            ipc_msg.err_code = (tri_data_sz == -1u ?
                    ERR_TRIANGULATE_GENERIC : ERR_INVALID_POLYNOMIAL);
            goto Lexit;
        }
        memcpy(&reply->surface_ID, &surface_ID, sizeof(MD5_digest));

//         write_log("SurfaceID: %u %u %u %u\n", surface_ID.md[0],
//                   surface_ID.md[1], surface_ID.md[2], surface_ID.md[3]);
        ipc_msg.shm_size = sizeof(SHM_Analysis_reply) + tri_data_sz;
        
        // remove client from the active list
        pthread_mutex_lock(&active_cl_mtx);
        active_clients.erase(host_addr.s_addr); 
        pthread_mutex_unlock(&active_cl_mtx);
        
    } else if(type == COMMENT) { // for PING request no special processing is
                                 // required => just send it back
        if(!handle_comment_request(&request_info, poly, host_addr)) {
            ipc_msg.err_code = ERR_REQUEST_PENDING;
            goto Lexit;
        }
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
    //pthread_cleanup_pop(do_cleanup); 
    return NULL;
}


void XTri_server::run()
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
        if(ipc_msg.m_type == ANALYSE || ipc_msg.m_type == TRIANGULATE ||
            ipc_msg.m_type == COMMENT || ipc_msg.m_type == PING ) {
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

// edit /etc/init.d/xtri_server
//install using: update-rc.d xtri_server defaults

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

    static XTri_server server_obj(inst);
    static Skeletonizer skeletonizer_obj;
//    QApplication q_app(argc, argv); 
    pxtri_server = &server_obj;
    pskeletonizer = &skeletonizer_obj;

    pxtri_server->setup();
    pxtri_server->run();
    return 0;
}

