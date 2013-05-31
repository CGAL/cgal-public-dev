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
// File          : xalci_client.C
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:04 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#include <signal.h>
#include <stdarg.h> 
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>

#include <errno.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/shm.h>
#include <pthread.h>
#include <semaphore.h>

#include <arpa/inet.h>
#include <set>
#include <string>
#include <math.h>
#include <memory.h>

#include "include/xtri_client.h"
#include "include/cgic.h"

static CGI_Client my_app;
static char static_poly[BUF_SIZE];

pthread_attr_t detach_attr;
sem_t shadow_sem;

bool CGI_Client::get_params(Request_type req, Triangulate_mode mode) {

    if(req == REQ_COMMENT) {
        if(cgiFormString("comment",static_poly,2001) != cgiFormSuccess)
            return false;
        return true;
    }

    cgiFormDouble("left", &in_shm_data.en_left, 0.05);
    cgiFormDouble("right", &in_shm_data.en_right, 0.05);
    cgiFormDouble("btm", &in_shm_data.en_bottom, 0.05);
    cgiFormDouble("top", &in_shm_data.en_top, 0.05);

    cgiFormDouble("below", &in_shm_data.z_below, -2);
    cgiFormDouble("above", &in_shm_data.z_above, 2);
    
    cgiFormInteger("sx", &in_shm_data.sx, 1);
    cgiFormInteger("sy", &in_shm_data.sy, 1);
    
    // in case of space subdivision we are passed the polynomial
    // equation and not an MD5 checksum
    if(req == REQ_ANALYSE ) {
        if(cgiFormString("surface", static_poly, BUF_SIZE) !=
                cgiFormSuccess)
            return false;
        
    } else if(req == REQ_TRIANGULATE) {
        char buf[256], *token;
        if(cgiFormString("surfaceid", buf, 256) != cgiFormSuccess)
            return false;
        token = strtok(buf, " ");
        int i=0;
        while(token != NULL && i < 4) {
            poly_hash.md[i++] = atoi(token);
            token = strtok(NULL, " ");
        }
    }
    return true;
}

void *message_thread(void *data)
{
    Thread_info *info = (Thread_info *)data;
    if(msgrcv(my_app.get_msg_queue_id(), info->pmsg, sizeof(IPC_Message)-4, 
            info->response, 0) == -1) {
        err_msg("msgrcv");
        err_exit();
    }
    sem_post(&shadow_sem); // unlock semaphore - response received
    return NULL;
}

bool CGI_Client::process(Request_type req, Triangulate_mode mode,
                         int server_id)
{

    server_id %= N_SERVER_INSTANCES;
    if(server_id < 0 || server_id >= N_SERVER_INSTANCES)
        server_id = 0;

    key_t key;
    key = ftok(KEY_FILENAME, 'm') + server_id + WEBXTI_UNIQUE_KEY;
    if((mq_id = msgget(key, IPC_CREAT|0666)) == -1) {
        err_msg("msgget");
        err_exit();
    }
    err_code = ERR_INVALID_DATA;
    //fprintf(stderr, "referrer: %s\n", cgiReferrer);

    if(req == REQ_COMMENT) {
        if(strstr(cgiReferrer, "/feedback.html") == NULL) {
            err_code = ERR_INVALID_REFERRER;
            return false;
        }
    } /*else if(strstr(cgiReferrer, "cgi-bin/xtri.cgi") == NULL) {

        fprintf(stderr, "invalid referrer: %s\n",
                cgiReferrer);
        // check validity code
        char temp[64];
        if(cgiFormString("x", temp, 32) != cgiFormSuccess
            || strcmp(temp, "Tg1872GHY67UT098oo") != 0) {
            err_code = ERR_INVALID_REFERRER;
            return false;
        }
    }*/
    //cgiHeaderCookieSetInteger("timestamp", (int)now, 86400*365, "/", "");
    int n_tries = 50;
    key = (key_t)(time(NULL) + PID);
    while(1) { // attempt to create a new key descriptor
        if((shm_id = shmget(key, MAX_SEG_SIZE,
                    IPC_CREAT|IPC_EXCL|0660))!= -1)
            break;
        if(errno != EEXIST) {
            err_msg("shmget");
            err_exit();
        }
        if(--n_tries <= 0) {
            fputs("shmget: failed to create a shared memory segment..\n",
                stderr);
            err_exit();
        }
        key++;
    }
    if((int&)(shm_addr = (char *)shmat(shm_id, NULL, 0)) == -1) {
        err_msg("shmat");
        err_exit();
    }
    
    IPC_Message ipc_msg;
    ipc_msg.shm_key = key;
    ipc_msg.shm_size = sizeof(SHM_Data);
    int resp;
    Message_type msg_type;
    
    SHM_Data *data = (SHM_Data *)shm_addr;
    // copy the input parameters
    memcpy(data, &in_shm_data, sizeof(SHM_Data));
   
    data->PID = PID;
    data->host_addr = inet_addr(cgiRemoteAddr);
    data->mode = mode;
    int *ascii_poly = (int *)(shm_addr + sizeof(SHM_Data));

    // check client data for sanity
    if(req == REQ_ANALYSE || req == REQ_TRIANGULATE) {
      
        bool failed = (data->z_below >= data->z_above ||
           fabs(data->z_below) > 50.0 || fabs(data->z_above) > 50.0 ||
           data->sx < 1 || data->sx > 50 || data->sy < 1 || data->sy > 50);

        if(failed)
            return false;
           
        if(mode == TRIANGULATE_DEFAULT) {
            failed = (data->en_left < 0.0 || data->en_left > 2.0 ||
                data->en_right < 0.0 || data->en_right > 2.0 ||
                data->en_bottom < 0.0 || data->en_bottom > 2.0 ||
                data->en_top < 0.0 || data->en_top > 2.0);
        } else if(mode == TRIANGULATE_ABS_BOUNDS) {
            failed = (fabs(data->en_left) > 50.0 || fabs(data->en_right) > 50.0
              || fabs(data->en_bottom) > 50.0 || fabs(data->en_top) > 50.0
              || data->en_left >= data->en_right
              || data->en_bottom >= data->en_top);
        } else  // bogus triangulation mode
            return false;

        if(failed) {
            fputs("wrong client data\n", stderr);
            return false;
        }
    }

    uint timeout = CLIENT_TIMEOUT;
    switch(req) {
    case REQ_ANALYSE:
        msg_type = ANALYSE;
        resp = (PID << 4)|ANALYSE_ACK;
        ipc_msg.shm_size += strlen(static_poly) + 1;
        strcpy((char *)ascii_poly, static_poly);
        break;
        
    case REQ_TRIANGULATE:
        ipc_msg.shm_size += sizeof(MD5_digest);
        memcpy((char *)ascii_poly, &poly_hash, sizeof(MD5_digest));
        msg_type = TRIANGULATE;
        resp = (PID << 4)|TRIANGULATE_ACK;
        break;
        
    case REQ_COMMENT:
        msg_type = COMMENT;
        resp = (PID << 4)|COMMENT_ACK;
/*        ipc_msg.shm_size += strlen(static_poly)+1;
        strcpy((char *)indices, static_poly);*/
        break;

    case REQ_PING:
        msg_type = PING;
        resp = (PID << 4)|PING_ACK;
        timeout = PING_TIMEOUT;
        break;
        
    default:
        return false;
    }
    
    ipc_msg.m_type = msg_type;
    if(msgsnd(mq_id, &ipc_msg, sizeof(IPC_Message)-4, 0) == -1) {
        err_msg("msgsnd");
        err_exit();
    }
    
    Thread_info info;
    pthread_t msg_thread_id;
    timespec ts;
    info.pmsg = &ipc_msg;
    info.response = resp;
        
    sem_init(&shadow_sem, 0, 0);
    pthread_attr_init(&detach_attr);
    pthread_attr_setdetachstate(&detach_attr, PTHREAD_CREATE_DETACHED);
    pthread_create(&msg_thread_id, &detach_attr, message_thread, &info);
        
    clock_gettime(CLOCK_REALTIME, &ts);
    ts.tv_sec += timeout; 
    if(sem_timedwait(&shadow_sem, &ts)!=0) 
        if(errno == ETIMEDOUT) {
            fputs("No connection to the server: cancelling receipt\n", stderr);
            pthread_cancel(msg_thread_id);  
            err_code = ERR_SERVER_TIMEOUT;
            return false;
        } else {
            err_msg("sem_timedwait");
            err_exit();
        }

    fprintf(stderr, "%d server response received.. forwarding\n", req);
        
    if(((unsigned *)shm_addr)[0] == 0xdeadbeef &&
            ((unsigned *)shm_addr)[1] == 0xdeadbeef) {
        err_code = ipc_msg.err_code;
        return false;
        
    } else {
        if(ipc_msg.shm_size < 0||ipc_msg.shm_size > MAX_SEG_SIZE) {
            return false;
        } 
        if(msg_type == ANALYSE || msg_type == TRIANGULATE) {

            SHM_Analysis_reply *reply = (SHM_Analysis_reply *)shm_addr;

            fprintf(cgiOut, "Content-type: application/octet-stream; charset=x-user-defined\r\n\r\n");
            
/*            fprintf(cgiOut, "Content-type: text/plain; charset=x-user-defined\r\n\r\n");*/
            
            fwrite(&reply->surface_ID, sizeof(MD5_digest), 1, cgiOut);
            // copy whatever other data you want to pass to the client
            fwrite(shm_addr + sizeof(SHM_Analysis_reply),
                    ipc_msg.shm_size - sizeof(SHM_Analysis_reply), 1, cgiOut);
        }
    }
    return true;
}

void CGI_Client::send_error_code() {

    //cgiHeaderContentType("image/png");
    fprintf(cgiOut, "Content-type: text/plain; charset=x-user-defined\r\n\r\n");

    // surface ID 0xfffff... indicates an error response
    MD5_digest err;
    memset(&err, 0xff, sizeof(MD5_digest));
    fwrite(&err, sizeof(MD5_digest), 1, cgiOut);
    fwrite(&err_code, sizeof(err_code), 1, cgiOut);
}

void CGI_Client::run()
{
    PID = getpid();
    Request_type req;
    Triangulate_mode mode = TRIANGULATE_DEFAULT;
    int sid = -1;

    //fputs("request received..\n", stderr);
    cgiFormInteger("id", (int *)&req, REQ_ANALYSE);
    cgiFormInteger("r", (int *)&mode, TRIANGULATE_DEFAULT);
    // get server id to communicate with
    cgiFormInteger("sid", (int *)&sid, -1);

    switch(req) {
    case REQ_ANALYSE:
    case REQ_TRIANGULATE:
        if(!get_params(req, mode) || !process(req, mode, sid)) {
            send_error_code();
        }
        break;

/*    case REQ_COMMENT:
        if(!get_params(req, mode) || !process(req, mode)) {
            send_error_code();
        }
        break;*/
        
    case REQ_PING: {
        cgiHeaderContentType("text/html");

        srand(time(NULL));
        if(sid == -1) {
            sid = rand() % N_SERVER_INSTANCES;
        }
            
        int start_sid = sid;
        while(1) {
            if(process(req, mode, sid)) {// check if this server is alive
                msg("error=%d&sid=%d", ERR_OK, sid); // save its ID
                break;
            }
            sid = (sid + 1) % N_SERVER_INSTANCES;
            if(sid == start_sid) { // all servers enumerated
                fprintf(stderr, "no response from servers\n");
                msg("error=%d", err_code);
                break;
            }
        }
        break;
    }    
    default:
        send_error_code();
        return;
    } 
}

int main(int argc, char *argv[]) 
{
    if(cgiInit(argc, argv)!=0)
        return -1;

    signal(SIGINT, sigint_handler);
    my_app.run();
    return 0;
}

std::string CGI_Client::parse_output(char *poly)
{
    std::string res(poly);
    while(1) {
        std::size_t idx = res.find("<");
        if(idx == std::string::npos)
            break;
        res.replace(idx, 1, "&lt;");
    }
    return res;
}

void sigint_handler(int)
{
    my_app.~CGI_Client();
    signal(SIGINT, SIG_DFL);
    raise(SIGINT);
}

CGI_Client::~CGI_Client()
{
    if(shm_id != -1) {
        shmdt(shm_addr);
        shmctl(shm_id, IPC_RMID, NULL);
    }
    sem_post(&shadow_sem);
    sem_destroy(&shadow_sem);
    pthread_attr_destroy(&detach_attr);
    cgiFreeResources();
}

void err_exit()
{
    exit(1);
}

void err_msg(const char *text)
{
    fprintf(stderr, "File: %s; line: %d\n", __FILE__, __LINE__);
    perror(text);
}

// CGI OUT formatted
void msg(const char *format, ...)
{
    va_list argptr;
    va_start(argptr, format);
    vfprintf(cgiOut, format, argptr);
    va_end(argptr);
}

// CGI ECHO
void echo(const char *text)
{
    fputs(text, cgiOut);
}
